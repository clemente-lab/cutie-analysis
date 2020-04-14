import pandas as pd
import numpy as np
import glob
np.random.seed(0)
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import click
from collections import defaultdict


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version='0.1')

# Required arguments
@click.option('-fv', '--fold_value', type=str,
              help='fold value for criterion for p value change')
@click.option('-s', '--statistic', type=str,
              help='string denoting type of analysis')
@click.option('-m', '--multi_corr', type=str,
              help='string denoting type of multiple corrections')
@click.option('-p', '--param', type=str,
              help='string denoting params used')
@click.option('-c', '--corr_compare', type=str,
              help='boolean denoting whether performing cooksd or not')
@click.option('-i', '--input_dir', type=click.Path(exists=True),
              help='input dir with .txt files of data')
@click.option('-o', '--output_dir', type=click.Path(exists=True),
              help='output dir to put config files')


def analyze_simulations_real(fold_value, statistic, multi_corr, param,
                             corr_compare, input_dir, output_dir):
    '''
    Script for analysis of real data by CUTIE
    '''

    def parse_log(f, cookd):
        lines = [l.strip() for l in f.readlines()]
        defaulted = False
        if cookd == 'True':
            for l in lines:
                #if "number of" in l:
                #    n_corr = int(l.split(' ')[-1])
                # el
                if "defaulted" in l:
                    defaulted = True
                elif "initial_corr" in l:
                    initial_corr = int(l.split(' ')[-1])
                elif "false correlations according to cookd" in l:
                    false_corr = int(l.split(' ')[-1])
                elif "true correlations according to cookd" in l:
                    true_corr = int(l.split(' ')[-1])
                elif "runtime" in l:
                    runtime = float(l.split(' ')[-1])
            rs_false = 0
            rs_true = 0

        else:
            # check if FDR correction defaulted
            for l in lines:
                # if "number of" in l:
                #    n_corr = int(l.split(' ')[-1])
                #el
                if "defaulted" in l:
                    defaulted = True
                elif "initial_corr" in l:
                    initial_corr = int(l.split(' ')[-1])
                elif "false correlations" in l:
                    false_corr = int(l.split(' ')[-1])
                elif "true correlations" in l:
                    true_corr = int(l.split(' ')[-1])
                elif "FP/TN1" in l:
                    rs_false = int(l.split(' ')[-1])
                elif "TP/FN1" in l:
                    rs_true = int(l.split(' ')[-1])
                elif "runtime" in l:
                    runtime = float(l.split(' ')[-1])

        # return n_corr, defaulted, initial_corr, false_corr, true_corr, rs_false, rs_true, runtime
        return defaulted, initial_corr, false_corr, true_corr, rs_false, rs_true, runtime


    headers = [
        'analysis_id',
        'parameter',
        'dataset',
        'statistic',
        'mc_used',
        'fold_value',
        'cooksd',
        # 'n_corr',
        'initial_corr',
        'true_corr(TP_FN)',
        'false_corr(FP_TN)',
        'rs_true_corr_TP_FN',
        'rs_false_corr_FP_TN',
        'true_frac',
        'false_frac',
        'rs_true_frac',
        'runtime']


    # populate df
    results_df = pd.DataFrame()

    mcs = multi_corr.split(',')
    fvs = fold_value.split(',')
    stats = statistic.split(',')
    cds = corr_compare.split(',')
    ds = ['lungc', 'lungpt', 'lungtx','who','hdac']
    params = param.split(',')
    for p in params:
        for mc in mcs:
            for fv in fvs:
                for s in stats:
                    for cd in cds:
                        for d in ds:
                            # nomc_10_pearson_True_lungpt
                            analysis_id = '_'.join([p, mc, fv, s, cd, d])
                            path = input_dir + analysis_id + '/'
                            files = sorted(glob.glob(path + '*.txt'))
                            # grab most recent log file
                            try:
                                rel_logfile = files[-1]
                                with open(rel_logfile, 'r') as f:
                                    # n_corr,
                                    defaulted, initial_corr, false_corr, \
                                        true_corr, rs_false, rs_true, runtime = parse_log(f,cd)

                                    true_frac = true_corr / initial_corr
                                    false_frac = false_corr / initial_corr
                                    rs_true_frac = rs_true / initial_corr

                                    new_row = pd.DataFrame([[analysis_id, p, d, s,
                                                            mc, fv, cd, # n_corr,
                                                            initial_corr, true_corr,
                                                            false_corr, rs_true,
                                                            rs_false, true_frac,
                                                            false_frac, rs_true_frac,
                                                            runtime]],
                                                            columns=headers)

                                    results_df = results_df.append(new_row)
                            except:
                                print(analysis_id)
                                print('Failed parsing')
                                if cd == 'True':
                                    if s == 'pearson':
                                        print(analysis_id)
                                else:
                                    print(analysis_id)

    results_df.to_csv(output_dir + 'real_results_df.txt', sep='\t')


    datasets = ['LungCancer', 'LungTranscriptomics', 'Gene Expression', 'WHO']

    col_to_id = {
        'LungCancer': 'lungc',
        'LungTranscriptomics': 'lungtx',
        'Gene Expression': 'hdac',
        'WHO': 'who'
    }

    statistics = ['pearson', 'spearman', 'kendall']


    # populate indices and ids for the dataframe and barplot
    for p in params:
        for fv in fvs:
            for mc in mcs:
                # subset to get relevant dataframe
                df = results_df[results_df['parameter'] == p]
                df = df[df['fold_value'] == fv]
                df = df[df['mc_used'] == mc]

                # create analysis id, e.g. p_fdr_1
                analysis_id = '_'.join([p, mc, fv])

                # create figure
                fig, axarr = plt.subplots(nrows=1, ncols=len(datasets),
                                          figsize=(10,6))

                # Custom x axis
                plt.xlabel("Dataset")

                # generate subplot x ticks
                x_stats = ['r', '\u03C1', '\u03C4'] # latter are rho and tau

                # set labels
                labels = ['True Positive', 'reverse sign-True Positive',
                          'False Positive', 'False Negative', 'True Negative']

                # set colors
                colors = ['#66b3ff','#ADD8E6','#ff9999','#99ff99','#8064A2']

                # iterate over datasets
                first = True
                for d, name in enumerate(datasets):
                    # iterate over statistic
                    stat_to_vals = defaultdict(list)

                    # hold values for stacked barplot
                    stat_to_sizes = {}

                    for s, stat in enumerate(statistics):
                        # extend analysis id, e.g. p_fdr_1_spearman_False_hdac
                        for_analysis_id = '_'.join([analysis_id, stat, 'False', col_to_id[name]])
                        rev_analysis_id = '_'.join([analysis_id, 'r' + stat, 'False', col_to_id[name]])

                        # get two relevant entries of df
                        for_df = df[df['analysis_id'] == for_analysis_id]
                        rev_df = df[df['analysis_id'] == rev_analysis_id]

                        # labels = ['TP', 'rsTP', 'FP', 'FN', 'rsFN', 'TN']

                        # TP is blue FP is red FN is green TN is purple
                        # for rs case
                        # reverse sign but still true FP is non reverse sign
                        # grab proportions of TP, rsTP, etc.
                        P = for_df['initial_corr'].values[0]
                        N = rev_df['initial_corr'].values[0]

                        # total = for_df['n_corr'].values[0]
                        total = P + N

                        TP = for_df['true_frac'].values[0]
                        rsTP = for_df['rs_true_frac'].values[0]

                        FN = rev_df['true_frac'].values[0]
                        rsFN = rev_df['rs_true_frac'].values[0]

                        # sizes = [(TP - rsTP) * P, rsTP * P,(1-TP)*P, (FN - rsFN) * N, rsFN * N, (1-FN)*N]
                        sizes = [(TP - rsTP) * P, rsTP * P,(1-TP)*P, FN * N, (1-FN)*N] / total
                        stat_to_sizes[stat] = sizes

                    # create df
                    raw_data = defaultdict(list)
                    for j, label in enumerate(labels):
                        for s in statistics:
                            raw_data[label].append(stat_to_sizes[s][j])

                    raw_df = pd.DataFrame(raw_data)

                    # set number of bars (# of statistics)
                    r = range(len(statistics))

                    # define subplot
                    ax = plt.subplot(1, len(datasets), d+1)

                    # ensure white background per plot
                    sns.set_style('white')

                    # build bottom bar stack
                    # fig = plt.figure(figsize=(8,4))
                    complete = np.zeros(len(statistics))
                    for k, label in enumerate(labels):
                        # create bars
                        plt.bar(r, raw_df[label], bottom = complete, color=colors[k],
                                edgecolor='white', width=0.85, label=label)
                        complete = np.add(complete, raw_df[label])

                    # subplot x ticks
                    plt.xticks(r, x_stats)

                    # remove axes
                    ax.spines['right'].set_visible(False)
                    ax.spines['top'].set_visible(False)

                    # remove y ticks if not first plot
                    if not first:
                        plt.yticks([])
                        ax.spines['left'].set_visible(False)

                    first = False

                    # dataset x label
                    plt.xlabel(name)


                # Add a legend: not useful when rs-TP is negligible
                # plt.legend(labels, loc='right', bbox_to_anchor=(1,1), ncol=1)

                # ensures legend is not cropped
                # plt.tight_layout()

                # save and close figure
                fig.savefig(output_dir + 'barplots_' + p + '_' + mc + '_' + fv + '.pdf')
                plt.close(fig)

    # generate figure 2
    # specific strings
    # p_fdr_1_kendall_False_lungc
    # p_fdr_3_pearson_False_hdac
    # p_fdr_1_spearman_False_lungtx
    # p_fdr_1_spearman_False_who

    fig_datasets = ['LC', 'LT', 'GE', 'WHO']

    ds_to_analyses = {
        'LC': ['p_fdr_1_kendall_False_lungc','p_fdr_1_rkendall_False_lungc'],
        'LT': ['p_fdr_1_spearman_False_lungtx','p_fdr_1_rspearman_False_lungtx'],
        'GE': ['p_fdr_3_pearson_False_hdac','p_fdr_3_rpearson_False_hdac'],
        'WHO': ['p_fdr_1_spearman_False_who','p_fdr_1_rspearman_False_who'],

    }

    analyses = []
    for ds in ds_to_analyses:
        analyses.extend(ds_to_analyses[ds])

    df = results_df[results_df['analysis_id'].isin(analyses)]

    # '\u03C1' is rho, '\u03C4' is tau
    x_labels = ['LC (\u03C4)',# fv = 1',
                'LT (\u03C1)',# fv = 1',
                'GE (r)',# fv = 3',
                'WHO (\u03C1)'#, fv = 1']

    # iterate over datasets
    ds_to_vals = defaultdict(list)

    # hold values for stacked barplot
    ds_to_sizes = {}

    for d, ds in enumerate(fig_datasets):
        # extend analysis id, e.g. p_fdr_1_spearman_False_hdac
        for_analysis_id, rev_analysis_id = ds_to_analyses[ds]

        # get two relevant entries of df
        for_df = df[df['analysis_id'] == for_analysis_id]
        rev_df = df[df['analysis_id'] == rev_analysis_id]

        # total = for_df['n_corr'].values[0]
        P = for_df['initial_corr'].values[0]
        N = rev_df['initial_corr'].values[0]
        total = P + N

        # grab fractions
        TP = for_df['true_frac'].values[0]
        rsTP = for_df['rs_true_frac'].values[0]

        FN = rev_df['true_frac'].values[0]
        rsFN = rev_df['rs_true_frac'].values[0]

        sizes = [(TP - rsTP) * P, rsTP * P,(1-TP)*P, FN * N, (1-FN)*N] / total
        ds_to_sizes[ds] = sizes

    # create df
    raw_data = defaultdict(list)
    for j, label in enumerate(labels):
        for ds in fig_datasets:
            raw_data[label].append(ds_to_sizes[ds][j])

    raw_df = pd.DataFrame(raw_data)

    # set number of bars (# of datasets)
    r = range(len(fig_datasets))

    # ensure white background per plot
    sns.set_style('white')

    # build bottom bar stack
    fig = plt.figure(figsize=(6,4))
    complete = np.zeros(len(fig_datasets))
    for k, label in enumerate(labels):
        # create bars
        plt.bar(r, raw_df[label], bottom = complete, color=colors[k],
                edgecolor='white', width=0.75, label=label)
        complete = np.add(complete, raw_df[label])

    # subplot x ticks
    plt.xticks(r, x_labels)

    # Custom x axis
    plt.xlabel("Dataset")

    # remove axes
    sns.despine()

    plt.tight_layout()
    fig.savefig(output_dir + 'Figure2_raw.pdf')
    plt.close(fig)


if __name__ == "__main__":
    analyze_simulations_real()




