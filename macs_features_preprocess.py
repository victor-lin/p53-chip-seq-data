import argparse
import pandas as pd


def output_features(data_file, minsamples, maxrep, output_fp):
    df_all = pd.read_table(data_file)
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
    df_subset.drop(['id', 'sample_count_distinct'], axis=1, inplace=True)
    df_subset.to_csv(output_fp, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_file', required=True,
                        help='tab-delimited data file')
    parser.add_argument('--minsamples', required=False, type=int, default=1,
                        help='minimum number of samples for binding intervals to be included (does not apply to non-binding intervals)')
    parser.add_argument('--maxrep', required=False, type=float, default=1.0,
                        help='maximum proportion of repeats for binding intervals to be included (does not apply to non-binding intervals)')
    parser.add_argument('-o', required=True,
                        help='Output filepath')
    args = parser.parse_args()
    output_features(args.data_file, args.minsamples, args.maxrep, args.o)
