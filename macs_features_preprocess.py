import argparse
import pandas as pd
import numpy as np

np.random.seed(0)


def output_features(data_file, minsamples, maxrep, test_size, train_out, test_out):
    df_all = pd.read_table(data_file)

    # subset data based on input conditions
    rows_binding_drop_n_samples = (df_all['sample_count_distinct'] > 0) & (df_all['sample_count_distinct'] < minsamples)
    rows_binding_drop_repeat_prop = (df_all['sample_count_distinct'] > 0) & (df_all['repeat_proportion'] > maxrep)
    rows_binding_drop = rows_binding_drop_n_samples | rows_binding_drop_repeat_prop

    ids_drop = set(df_all.loc[rows_binding_drop, 'id'].values)
    rows_all_drop = df_all['id'].isin(ids_drop)

    log = ('{} intervals ({} binding, {} non-binding) will be dropped out of '
           '{}.'.format(sum(rows_all_drop),
                        len(ids_drop),
                        sum(rows_all_drop) - len(ids_drop),
                        len(df_all)))
    print log

    df_subset = df_all[~rows_all_drop].copy()
    df_subset['max_MACS_score'] = (df_subset['max_MACS_score'] > 0).astype(int)
    df_subset.rename(columns={'max_MACS_score': 'binding'}, inplace=True)

    # generate train/test split
    ids_remaining = list(set(df_subset['id'].values))
    num_test = int(len(ids_remaining) * test_size)
    ids_test = np.random.choice(ids_remaining, size=num_test, replace=False)
    df_test = df_subset[df_subset['id'].isin(ids_test)].copy()
    df_test.drop(['id', 'sample_count_distinct'], axis=1, inplace=True)
    df_test.to_csv(test_out, sep='\t', index=False)

    df_train = df_subset[~df_subset['id'].isin(ids_test)].copy()
    df_train.drop(['id', 'sample_count_distinct'], axis=1, inplace=True)
    df_train.to_csv(train_out, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_file', required=True,
                        help='tab-delimited data file')
    parser.add_argument('--minsamples', required=False, type=int, default=1,
                        help='minimum number of samples for binding intervals to be included (does not apply to non-binding intervals)')
    parser.add_argument('--maxrep', required=False, type=float, default=1.0,
                        help='maximum proportion of repeats for binding intervals to be included (does not apply to non-binding intervals)')
    parser.add_argument('--test_size', required=True, type=float,
                        help='Proportion of data to be used for testing set (0,1)')
    parser.add_argument('--train_out', required=True,
                        help='output file for training data (tab-delimited)')
    parser.add_argument('--test_out', required=True,
                        help='output file for testing data (tab-delimited)')
    args = parser.parse_args()
    output_features(args.data_file, args.minsamples, args.maxrep,
                    args.test_size, args.train_out, args.test_out)
