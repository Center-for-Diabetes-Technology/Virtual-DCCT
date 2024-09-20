# Author: Ben Lobo (lobo@virginia.edu)

import itertools as itr
import json
import numpy as np
from operator import itemgetter
import pandas as pd
from scipy.interpolate import CubicSpline


def load_rtss_set(path):
    """
    Loads the RTSS set from the path provided.
    """
    json_rtss_set = load(path)
    
    rtss_set = {'Gamma': json_rtss_set['Gamma'],
                'Tau': json_rtss_set['Tau'],
                'RTSSs': {}}
    for rtss_idx, json_rtss in json_rtss_set['RTSSs'].items():
        rtss = {'tss': pd.read_json(json_rtss['tss'], 
                                    orient='split', 
                                    typ='series'),   
                'rs_tss': pd.read_json(json_rtss['rs_tss'], 
                                       orient='split', 
                                       typ='series'),
                'UIDx': json_rtss['UIDx'],
                'UIDy': json_rtss['UIDy'], 
                'SID_SDTx': (json_rtss['SID_SDTx'][0], 
                             pd.to_datetime(json_rtss['SID_SDTx'][1])),
                'SID_SDTy': (json_rtss['SID_SDTy'][0], 
                             pd.to_datetime(json_rtss['SID_SDTy'][1])),                
                'score': json_rtss['score']
               }
        rtss['notnull'] = rtss['tss'].notnull()
           
        rtss_set['RTSSs'][int(rtss_idx)] = rtss

    return rtss_set


def load(path):
    """
    Wrapper for the json.load() function that
    includes a context manager.
    """
    with open(path, 'r') as read_file:
        data = json.load(read_file)

    return data


def generate_CGM_data(PID, PID_rng, smbg_data_path, db_path, output_path):
    """
    Takes the dataframe smbgs of daily SMBG profile data and
    generates a daily CGM profile for each daily SMBG profile.
    Writes the CGM data to output_path.
    """
    smbgs = pd.read_csv(f'{smbg_data_path}/pid{PID}_smbg_profile_data.csv',
                        parse_dates=['Date', 'Date_Time']) \
                .set_index('Date')
    dates = sorted(smbgs.index.unique())

    data = {'smbg': [],
            'smbg_uncensored': [],
            'motif': [],
            'indices': [],
            'dp': []}
    with sqlite3.connect(db_path) as conn:
        c = conn.cursor()

        for day_index, date in enumerate(dates):
            # print(day_index, date)
            smbg_profile = smbgs.loc[date] \
                                .set_index('Time_Index') \
                                ['SMBG']
            smbg_uncensored_profile = smbgs.loc[date] \
                                            .set_index('Time_Index') \
                                            ['SMBG_Uncensored']
            data['smbg'].append(smbg_time_series(day_index, smbg_profile))
            data['smbg_uncensored'].append(smbg_time_series(day_index, smbg_uncensored_profile))
            
            motif_index, dp_index, dp = select_CGM_profile(c, 
                                                           smbg_profile, 
                                                           motifs,
                                                           motif_dp_indices,
                                                           PID_rng)            
            
            data['motif'].append(motifs[motif_index]['tss'])
            data['indices'].append((motif_index, dp_index[0], dp_index[1]))
            data['dp'].append(dp)

    pd.DataFrame(data['indices'], 
                 columns=['Motif_Index', 'DP_SID', 'DP_Start_DT']) \
        .to_csv(f'{output_path}/pid{PID}_motif_dp_indices.csv', index=False)
                
    smbg_ts, smbg_uncensored_ts, motif_ts, dp_ts = data_to_time_series(data, dates)
    data = pd.concat([smbg_ts, smbg_uncensored_ts, motif_ts, dp_ts], axis=1) \
            .rename(columns={0: 'SMBG',
                             1: 'SMBG_Uncensored',
                             2: 'Motif',
                             3: 'Daily_Profile'})
    data.to_csv(f'{output_path}/pid{PID}_smbg_motif_dp.csv', 
                index=True)


def smbg_time_series(day_index, smbg_profile):
    """
    Creates a single day time series using the SMBG data.
    """    
    return pd.Series(smbg_profile.values, 
                     index=[ti + (day_index * 288) for ti in smbg_profile.index])


def select_CGM_profile(c, smbg_profile, motifs, motif_dp_indices, rng):
    """
    Randomly selects 5 daily CGM profiles for each motif in
    motif_indices, and selects the best match
    """
    motif_indices = matching_motifs(smbg_profile, motifs)

    matches = []
    for motif_index in motif_indices:
        dp_indices = sorted(motif_dp_indices[motif_index])
        rng.shuffle(dp_indices)

        num_dps = 0
        for SID, start_dt in dp_indices:
            dp = daily_CGM_profile(c, SID, start_dt)

            if dp.notnull().sum() > 230:  # 230 = 80% of data
                num_dps += 1
                matches.append({'ls_score': score(dp, smbg_profile),
                                'motif_index': motif_index,
                                'dp_index': (SID, start_dt),
                                'dp': dp})
            if num_dps >= 5:
                break

    # Sort the matches using 'ls_score'
    matches = sorted(matches, key=lambda d: d['ls_score'])
    
    return matches[0]['motif_index'], \
            matches[0]['dp_index'], \
            matches[0]['dp']


def matching_motifs(smbg_profile, motifs, threshold_pct=5):
    """
    Returns a list of motif indices (could be just one index)
    of the motifs which best match the SMBG profile passed
    in.
    """
    scores = []
    for motif_idx in range(0, 483):
        scores.append(score(motifs[motif_idx]['tss'], smbg_profile))
    scores = pd.Series(scores).sort_values(ascending=True)
    
    threshold = scores.iloc[0] * (1 + (threshold_pct / 100))
    
    return list(scores.loc[scores < threshold].index)


def score(ts, smbg_profile):
    """
    Returns the score (least squares estimate) between the time 
    series and the SMBG profile.
    A smaller score indicates that the time series matches the SMBG
    profile better.
    """
    # Check for any nan's in the ts which overlap the SMBG positions.  
    # If there are, interpolate the ts before proceeding.
    if len(set(ts.index[ts.isnull()]).intersection(smbg_profile.index)) > 0:
        ts = interpolate_missing_values(ts, 24)    
    
    lse = 0
    for time_index, smbg_val in smbg_profile.items():
        lse += (smbg_val - ts.loc[time_index]) ** 2
        
    return lse


def interpolate_missing_values(tss, tss_length_hrs):
    """
    Assumes the TSS contains data of a single patient of time 
    length tss_length_hrs.
    Interpolates missing values
    1. Using cubic splines if there are two data points at 
        the start and end of the missing data interval,
    2. Using linear interpolation (through cubic splines) 
        if there is only one data point available at either
        the start or end of the missing data interval, and 
    3. Using forward or backward fill as the interpolation 
        method for the interval if the missing data interval
        contains either the first or last data point in the 
        time series.
    """
    missing_intervals = interval_starts_and_length(tss.loc[tss.isnull()].index)
    if missing_intervals:
        if missing_intervals[0][0] == 0:
            missing_intervals.pop(0)
        if missing_intervals and \
            (missing_intervals[-1][0] + missing_intervals[-1][1] == (tss_length_hrs * 12)):
            missing_intervals.pop(-1)
            
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CubicSpline.html
        for start_tb, length in missing_intervals:
            # Try to use two points on either end of the missing 
            # data interval first, and if that is not possible
            # then only use one point on either end
            try:
                y = [tss.loc[start_tb - 2], tss.loc[start_tb - 1], 
                     tss.loc[start_tb + length], tss.loc[start_tb + length + 1]]
                
                if pd.isnull(y).any():
                    raise KeyError

                cs = CubicSpline(x=[start_tb - 2, start_tb - 1, 
                                    start_tb + length, start_tb + length + 1],
                                 y=y)
            except KeyError as e:
                cs = CubicSpline(x=[start_tb - 1, 
                                    start_tb + length],
                                 y=[tss.loc[start_tb - 1], 
                                    tss.loc[start_tb + length]])

            tbs_to_interpolate = range(start_tb, start_tb + length)
            tss.loc[tbs_to_interpolate] = sensor_bounding(cs(tbs_to_interpolate))

        return tss.fillna(method='ffill').fillna(method='backfill')
    else:
        return tss
    
    
def interval_starts_and_length(indices_list):
    """
    https://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
    https://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list
    """
    ranges = []
    for k, g in itr.groupby(enumerate(indices_list), lambda ix: ix[0] - ix[1]):
        ranges.append(list(map(itemgetter(1), g)))
        
    return [(ti_range[0], len(ti_range)) for ti_range in ranges]   


def sensor_bounding(BG_vals):
    """
    Makes sure that the list of BG values passed in
    conforms to the sensor bounds of 40 g/dL and 
    400 g/dL.
    """
    bounded_BG_vals = []
    for BG in BG_vals:
        if BG < 40:
            bounded_BG_vals.append(39)
        elif BG > 400:
            bounded_BG_vals.append(401)
        else:
            bounded_BG_vals.append(BG)
            
    return bounded_BG_vals


def daily_CGM_profile(c, SID, start_dt):
    """
    """  
    query = """
        SELECT dt_index, bg_mg_per_dL FROM blood_glucose
        WHERE SID == :SID
        AND datetime(dt_index) >= :start_dt
        AND datetime(dt_index) < datetime(:start_dt, '+' || :tss_length_hrs || ' hours')
        """
    dp_list = c.execute(query, {'SID': SID,
                                'start_dt': str(start_dt),
                                'tss_length_hrs': 24}) \
                .fetchall()

    dp = pd.DataFrame(dp_list, columns=['DT_Index', 'Value']) \
            .set_index('DT_Index') \
            ['Value']
    dp.index = pd.to_datetime(dp.index)

    return complete_tss(dp, start_dt, 24)


def complete_tss(tss, start_dt, tss_length_hrs):
    """
    Returns a pandas series with tss_length_hrs * 12 entries, 
    and the entries that were missing in the TSS passed
    in are set to NaN.
    Note that the index of the series is just the integers
    in range(tss_length_hrs * 12).
    """
    expected_n = tss_length_hrs * 12
    actual_n = tss.shape[0]

    tss = tss.sort_index()
    if actual_n < expected_n:
        expected_index = pd.date_range(start_dt, periods=12 * tss_length_hrs, freq='5T')
        return tss.reindex(expected_index) \
                    .reset_index(drop=True)
    else:
        return tss.reset_index(drop=True)
    
    
def data_to_time_series(data, dates):
    """
    """
    smbg_ts = pd.concat(data['smbg'], axis=0)
    smbg_uncensored_ts = pd.concat(data['smbg_uncensored'], axis=0)
    motif_ts = pd.concat(data['motif'], axis=0).reset_index(drop=True)
    dp_ts = pd.concat(data['dp'], axis=0).reset_index(drop=True)  
    
    dt_index = pd.date_range(dates[0], 
                             dates[-1] + pd.DateOffset(days=1), 
                             freq='5T', 
                             inclusive='left', 
                             name='Date_Time')
    return pd.Series(smbg_ts.reindex(range(motif_ts.shape[0])).values, index=dt_index), \
            pd.Series(smbg_uncensored_ts.reindex(range(motif_ts.shape[0])).values, index=dt_index), \
            pd.Series(motif_ts.values, index=dt_index), \
            pd.Series(dp_ts.values, index=dt_index)
