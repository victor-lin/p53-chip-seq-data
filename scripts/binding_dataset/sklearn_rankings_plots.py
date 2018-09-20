import os
import pandas as pd

from sklearn_plot_helper_functions import (
    get_svc_coefficients,
    get_svc_rfe_rankings,
    get_all_weights_plot,
    get_motif_weights_plot,
)

results_dir = '../../results'
figures_dir = os.path.join(results_dir, 'figures')
if not os.path.exists(figures_dir):
    print 'making directory {}'.format(figures_dir)
    os.makedirs(figures_dir)

# load non-rep data
fp_train_nonrep = '../../results/datafiles/peaks_merged_features_train__minsamples_2__rep_max0.1__nonbinding_20.txt'
fp_test_nonrep = '../../results/datafiles/peaks_merged_features_test__minsamples_2__rep_max0.1__nonbinding_20.txt'
df_train_nonrep = pd.read_table(fp_train_nonrep)
df_test_nonrep = pd.read_table(fp_test_nonrep)

# load rep data
fp_train_rep = '../../results/datafiles/peaks_merged_features_train__minsamples_2__rep_min0.9__nonbinding_20.txt'
fp_test_rep = '../../results/datafiles/peaks_merged_features_test__minsamples_2__rep_min0.9__nonbinding_20.txt'
df_train_rep = pd.read_table(fp_train_rep)
df_test_rep = pd.read_table(fp_test_rep)


def plot_svc_nonrep():
    weights = get_svc_coefficients(df_train_nonrep, df_test_nonrep)

    plt = get_motif_weights_plot(weights, df_train_nonrep.columns[1:], 0.7,
                                 'SVC motif feature weights (nonrepetitive)')
    save_path = os.path.join(figures_dir, 'motif_weights_svc_nonrep.png')
    plt.savefig(save_path, bbox_inches='tight')
    plt.clf()

    plt = get_all_weights_plot(weights, df_train_nonrep.columns[1:], 0.7,
                               'SVC most important features (nonrepetitive)')
    save_path = os.path.join(figures_dir, 'all_weights_svc_nonrep.png')
    plt.savefig(save_path, bbox_inches='tight')
    plt.clf()


def plot_svc_rep():
    weights = get_svc_coefficients(df_train_rep, df_test_rep)

    plt = get_motif_weights_plot(weights, df_train_rep.columns[1:], 0.7,
                                 'SVC motif feature weights (repetitive)')
    save_path = os.path.join(figures_dir, 'motif_weights_svc_rep.png')
    plt.savefig(save_path, bbox_inches='tight')
    plt.clf()

    plt = get_all_weights_plot(weights, df_train_rep.columns[1:], 0.7,
                               'SVC most important features (repetitive)')
    save_path = os.path.join(figures_dir, 'all_weights_svc_rep.png')
    plt.savefig(save_path, bbox_inches='tight')
    plt.clf()


def save_svc_rfe_rankings():
    # nonrep
    ranking_values = get_svc_rfe_rankings(df_train_nonrep, df_test_nonrep)
    feature_rankings = zip(ranking_values, df_train_nonrep.columns[1:])
    df = pd.DataFrame(zip(*sorted(feature_rankings)))
    df = df.transpose()
    df.columns = ['rank', 'feature_name']
    save_path = os.path.join(results_dir, 'svc_rfe_rankings_nonrep.txt')
    df.to_csv(save_path, index=False, sep='\t')

    # rep
    ranking_values = get_svc_rfe_rankings(df_train_rep, df_test_rep)
    feature_rankings = zip(ranking_values, df_train_rep.columns[1:])
    df = pd.DataFrame(zip(*sorted(feature_rankings)))
    df = df.transpose()
    df.columns = ['rank', 'feature_name']
    save_path = os.path.join(results_dir, 'svc_rfe_rankings_rep.txt')
    df.to_csv(save_path, index=False, sep='\t')


if __name__ == '__main__':
    plot_svc_nonrep()
    plot_svc_rep()
