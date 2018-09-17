import argparse
import pandas as pd
import numpy as np

from sklearn.preprocessing import MinMaxScaler

np.random.seed(0)


def output_features(data_file, test_size, train_out, test_out):
    df = pd.read_table(data_file)

    # set class label
    df['max_MACS_score'] = (df['max_MACS_score'] > 0).astype(int)
    df.rename(columns={'max_MACS_score': 'binding'}, inplace=True)

    # generate train/test split
    ids = list(set(df['id'].values))
    num_test = int(len(ids) * test_size)
    ids_test = np.random.choice(ids, size=num_test, replace=False)

    cols_X = list(df.columns)
    for col in ['binding', 'id', 'sample_count_distinct']:
        cols_X.remove(col)
    scaler = MinMaxScaler()

    df_test = df[df['id'].isin(ids_test)].copy()
    df_test.drop(['id', 'sample_count_distinct'], axis=1, inplace=True)
    df_test[cols_X] = scaler.fit_transform(df_test[cols_X])
    df_test.to_csv(test_out, sep='\t', index=False)

    df_train = df[~df['id'].isin(ids_test)].copy()
    df_train.drop(['id', 'sample_count_distinct'], axis=1, inplace=True)
    df_train[cols_X] = scaler.fit_transform(df_train[cols_X])
    df_train.to_csv(train_out, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_file', required=True,
                        help='tab-delimited data file')
    parser.add_argument('--test_size', required=True, type=float,
                        help='Proportion of data to be used for testing set (0,1)')
    parser.add_argument('--train_out', required=True,
                        help='output file for training data (tab-delimited)')
    parser.add_argument('--test_out', required=True,
                        help='output file for testing data (tab-delimited)')
    args = parser.parse_args()
    output_features(args.data_file, args.test_size, args.train_out, args.test_out)
