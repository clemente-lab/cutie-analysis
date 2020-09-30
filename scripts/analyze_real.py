import pandas as pd
import numpy as np
import glob
import os
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
@click.option('-d', '--datasets', type=str,
              help='string denoting datasets used')
@click.option('-c', '--corr_compare', type=str,
              help='boolean denoting whether performing cooksd or not')
@click.option('-i', '--input_dir', type=click.Path(exists=True),
              help='input dir with .txt files of data')
@click.option('-o', '--output_dir', type=click.Path(exists=True),
              help='output dir to put config files')


def analyze_real(fold_value, statistic, multi_corr, param, datasets,
                             corr_compare, input_dir, output_dir):
    '''
    Script for analysis of real data by CUTIE
    '''

    def parse_log(f, cookd):
        lines = [l.strip() for l in f.readlines()]
        defaulted = False
        if cookd == 'True':
            for l in lines:
                if "number of correlations" in l:
                    n_corr = int(l.split(' ')[-1])
                elif "is 0.05" in l:
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
                if "number of correlations" in l:
                    n_corr = int(l.split(' ')[-1])
                elif "is 0.05" in l:
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

        return n_corr, defaulted, initial_corr, false_corr, true_corr, rs_false, rs_true, runtime




    # lungc,lungtx,hdac,who
    # datasets = ['LungCancer', 'LungTranscriptomics', 'Gene Expression', 'WHO']
    # datasets = ['Microbiome', 'Gene Expression', 'World Health Statistics']

    # ds = ['lungc','hdac','who']
    # happens later
    # ds = datasets.split(',')

    # parallel_ds = ['livermfull', ...]
    # if d in []

    # '\u03C1' is rho, '\u03C4' is tau
    x_labels = ['LC (\u03C4)',# fv = 1',
                # 'LT (\u03C4)',# fv = 1',
                'GE (r)',# fv = 3',
                'WHO (\u03C1)']#, fv = 1']

    # more verbose labels
    x_labels = ['Lung\nCancer\n(\u03C4)',
            # 'Lung\nTranscriptomics\n(\u03C4)',
            'Gene\nExpression\n(r)',
            'WHO\n(\u03C1)']

    # temporary for testing
    x_labels = datasets.split(',')

    ds_to_titles = {
        'lungc': 'Lung Cancer (\u03C4)', # Microbiome
        'lungtx': 'Lung Transcriptomics (\u03C4)',
        'hdac': 'Gene Expression (r)',
        'baseball':'baseball',
        'micro':'micro',
        'microfull':'microfull',
        'who': 'WHO (\u03C1)', # 'World Health Statistics'
        'airplane': 'airplane',
        'airplanefull': 'airplanefull',
        'hdacfull': 'hdacfull',
        'lungpt': 'Lung Pneumotyping',
        'ad0': 'ad0',
        'ad1': 'ad1',
        'oom': 'oom',
        'roc': 'roc',
        'mennonites': 'mennonites',
        'covid': 'covid',
        'airplane': 'airplane',
        'ici':'ici',
        'liverf':'liverf',
        'liverm':'liverm',
        'hgoral': 'hgoral',
        'spatial':'spatial',
        'crc':'crc',
        'ibd': 'ibd',
        'cell':'cell',
        'nc':'nc',
        'plos':'plos',
        'ca':'ca',
        'statin':'statin',
        'livermfull':'livermfull',
        'liverffull':'liverffull'
    }

    statistics = statistic.split(',')

    for s in statistics:
        if s[0] == 'r':
            statistics.remove(s)

    # want statistics to be just forward
    # ['pearson', 'spearman', 'kendall', 'mine']

    headers = [
        'analysis_id',
        'parameter',
        'dataset',
        'statistic',
        'mc_used',
        'fold_value',
        'cooksd',
        'n_corr',
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
    datasets = datasets.split(',')
    params = param.split(',')

    if not os.path.exists(output_dir + 'real_results_df.txt'):
        for p in params:
            for mc in mcs:
                for fv in fvs:
                    for s in stats:
                        for cd in cds:
                            for d in datasets:
                                # construct identifier string
                                # nomc_10_pearson_True_lungpt
                                analysis_id = '_'.join([p, mc, fv, s, cd, d])

                                # for single jobs
                                try:
                                    path = input_dir + analysis_id + '/'
                                    files = sorted(glob.glob(path + '*.txt'))

                                    # grab most recent log file
                                    rel_logfile = files[-1]
                                    with open(rel_logfile, 'r') as f:
                                        n_corr, defaulted, initial_corr, false_corr, \
                                            true_corr, rs_false, rs_true, runtime = parse_log(f,cd)

                                        true_frac = true_corr / initial_corr
                                        false_frac = false_corr / initial_corr
                                        rs_true_frac = rs_true / initial_corr

                                        new_row = pd.DataFrame([[analysis_id, p, d, s,
                                                                mc, fv, cd, n_corr,
                                                                initial_corr, true_corr,
                                                                false_corr, rs_true,
                                                                rs_false, true_frac,
                                                                false_frac, rs_true_frac,
                                                                runtime]],
                                                                columns=headers)

                                        results_df = results_df.append(new_row)
                                # for split jobs
                                except:
                                    try:
                                        path = input_dir + analysis_id
                                        # if analysis id is in []
                                        split_jobs = glob.glob(path + '*')
                                        # nomc_10_pearson_True_lungpt_0/, etc.

                                        dfs = []
                                        for j in split_jobs:
                                            try:
                                                df = pd.read_csv(j + '/data_processing/summary_df_resample_1.txt', sep='\t')

                                                # filter out 'unpaired' identical var pairs
                                                df = df[df['var1'] != df['var2']]

                                                dfs.append(df)
                                            except:
                                                print(j)

                                        # merge all dfs
                                        final_df = pd.concat(dfs, axis=0)

                                        # sort columns to remove duplicate var pairs
                                        # final_df = final_df.sort_values(by=['var1', 'var2'])
                                        # final_df = final_df.drop_duplicates(subset=['var1', 'var2'], keep='last')
                                        # instead of trying to figure out how to drop duplicates, we divide length by 2
                                        # to account for double counting (assumes statistic is symmetric)

                                        n_corr = len(final_df) / 2
                                        initial_df = final_df[final_df['class'].isin(['TP','FP','TN','FN'])]
                                        initial_corr = len(initial_df) / 2
                                        true_df = initial_df[initial_df['class'].isin(['TP','FN'])]
                                        false_df = initial_df[initial_df['class'].isin(['FP','TN'])]
                                        true_corr = len(true_df) / 2
                                        false_corr = len(initial_df[initial_df['class'].isin(['FP','TN'])]) / 2
                                        rs_true = len(true_df[true_df['reverse'] == 'Yes']) / 2
                                        rs_false = len(false_df[false_df['reverse'] == 'Yes']) / 2

                                        true_frac = true_corr / initial_corr
                                        false_frac = false_corr / initial_corr
                                        rs_true_frac = rs_true / initial_corr
                                        runtime = 'parallel'

                                        new_row = pd.DataFrame([[analysis_id, p, d, s,
                                                                mc, fv, cd, n_corr,
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

        results_df.to_csv(output_dir + 'real_results_df.txt', sep='\t', index=False)
    else:
        results_df = pd.read_csv(output_dir + 'real_results_df.txt', sep='\t')
        for col in ['fold_value']:
            results_df[col] = results_df[col].astype(str)

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

                # ensure white background per plot with ticks
                sns.set(font_scale=1.2)
                sns.set_style("ticks", {'font.family':'sans-serif','font.sans-serif':'Helvetica'})

                # create figure
                fig, axarr = plt.subplots(nrows=1, ncols=len(datasets),
                                          figsize=(3 * len(datasets),6))

                # Custom x axis
                plt.xlabel("Dataset")

                # generate subplot x ticks
                x_stats = ['r', '\u03C1', '\u03C4'] # latter are rho and tau

                # set labels
                labels = ['True Positive', 'reverse sign-True Positive',
                          'False Positive', 'False Negative', 'True Negative']

                # set colors # #228B22  #99ff99
                colors = ['#66b3ff','#ADD8E6','#ff9999','#228B22','#8064A2']

                # iterate over datasets
                first = True
                for d, ds in enumerate(datasets):
                    # iterate over statistic
                    stat_to_vals = defaultdict(list)

                    # hold values for stacked barplot
                    stat_to_sizes = {}

                    for s, stat in enumerate(statistics):
                        # extend analysis id, e.g. p_fdr_1_spearman_False_hdac
                        for_analysis_id = '_'.join([analysis_id, stat, 'False', ds])
                        rev_analysis_id = '_'.join([analysis_id, 'r' + stat, 'False', ds])

                        # get two relevant entries of df
                        for_df = df[df['analysis_id'] == for_analysis_id]
                        rev_df = df[df['analysis_id'] == rev_analysis_id]

                        # labels = ['TP', 'rsTP', 'FP', 'FN', 'rsFN', 'TN']

                        # grab proportions of TP, rsTP, etc.
                        total = for_df['n_corr'].values[0]

                        P = for_df['initial_corr'].values[0] / total
                        N = rev_df['initial_corr'].values[0] / total

                        TP = for_df['true_frac'].values[0]
                        rsTP = for_df['rs_true_frac'].values[0]

                        FN = rev_df['true_frac'].values[0]
                        # rsFN = rev_df['rs_true_frac'].values[0]

                        # sizes = [(TP - rsTP) * P, rsTP * P,(1-TP)*P, (FN - rsFN) * N, rsFN * N, (1-FN)*N]
                        sizes = [(TP - rsTP) * P, rsTP * P,(1-TP)*P, FN * N, (1-FN)*N]
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

                    # build bottom bar stack
                    # fig = plt.figure(figsize=(8,4))
                    complete = np.zeros(len(statistics))
                    for k, label in enumerate(labels):
                        # create bars
                        plt.bar(r, raw_df[label], bottom = complete, color=colors[k],
                                edgecolor='white', width=0.85, label=label, linewidth=0)
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
                    plt.xlabel(ds)


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

    # analyses = []
    #for ds in ds_to_analyses:
    #    analyses.extend(ds_to_analyses[ds])

    # df = results_df[results_df['analysis_id'].isin(analyses)]
    df = results_df

    # iterate over datasets
    ds_to_vals = defaultdict(list)

    # hold values for stacked barplot
    ds_to_sizes = {}

    for d, ds in enumerate(datasets):
        # extend analysis id, e.g. p_fdr_1_spearman_False_hdac
        # for_analysis_id, rev_analysis_id = ds_to_analyses[ds]
        for_analysis_id = '_'.join(['p','nomc','1','pearson','False',ds])
        rev_analysis_id = '_'.join(['p','nomc','1','rpearson','False',ds])

        # get two relevant entries of df
        for_df = df[df['analysis_id'] == for_analysis_id]
        rev_df = df[df['analysis_id'] == rev_analysis_id]

        # get N and P
        total = for_df['n_corr'].values[0]
        P = for_df['initial_corr'].values[0] / total
        N = rev_df['initial_corr'].values[0] / total

        # grab fractions
        TP = for_df['true_frac'].values[0]
        rsTP = for_df['rs_true_frac'].values[0]

        FN = rev_df['true_frac'].values[0]
        rsFN = rev_df['rs_true_frac'].values[0]

        sizes = [(TP - rsTP) * P, rsTP * P,(1-TP)*P, FN * N, (1-FN)*N]
        ds_to_sizes[ds] = sizes

    # create df
    raw_data = defaultdict(list)
    for j, label in enumerate(labels):
        for ds in datasets:
            raw_data[label].append(ds_to_sizes[ds][j])

    raw_df = pd.DataFrame(raw_data, index=datasets)
    raw_df.to_csv(output_dir + 'Fig2_rawdata_df.txt', sep='\t')

    # set number of bars (# of datasets)
    r = range(len(datasets))

    # ensure white background per plot
    sns.set(font_scale=0.9)
    sns.set_style("ticks", {'font.family':'sans-serif','font.sans-serif':'Helvetica'})

    # build bottom bar stack
    fig = plt.figure(figsize=(5,4))
    complete = np.zeros(len(datasets))
    for k, label in enumerate(labels):
        # create bars
        plt.bar(r, raw_df[label], bottom = complete, color=colors[k],
                edgecolor='white', width=0.75, label=label, linewidth=0)
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

    #mcs = multi_corr.split(',')
    #fvs = fold_value.split(',')
    #stats = statistic.split(',')
    #cds = corr_compare.split(',')
    #ds = datasets.split(',')
    #params = param.split(',')

    # condensed df
    raw_data = defaultdict(list)
    for p in params:
        for mc in mcs:
            for fv in fvs:
                for cd in cds:
                    for d in datasets:
                        for stat in stats:
                            # do for only forward stats
                            if stat[0] != 'r':
                                for_analysis_id = '_'.join([p,mc,str(fv),stat,cd,d])
                                rev_analysis_id = '_'.join([p,mc,str(fv),'r'+stat,cd,d])

                                for_df = results_df[results_df['analysis_id'] == for_analysis_id]
                                rev_df = results_df[results_df['analysis_id'] == rev_analysis_id]

                                raw_data['analysis_id'].append(for_analysis_id)
                                raw_data['TP'].append(for_df['true_corr(TP_FN)'].values[0])
                                raw_data['FP'].append(for_df['false_corr(FP_TN)'].values[0])

                                if cd == 'False':
                                    raw_data['rsTP'].append(for_df['rs_true_corr_TP_FN'].values[0])
                                    raw_data['FN'].append(rev_df['true_corr(TP_FN)'].values[0])
                                    raw_data['TN'].append(rev_df['false_corr(FP_TN)'].values[0])
                                else:
                                    raw_data['rsTP'].append(0)
                                    raw_data['FN'].append(0)
                                    raw_data['TN'].append(rev_df['true_corr(TP_FN)'].values[0] + rev_df['false_corr(FP_TN)'].values[0])


    raw_df = pd.DataFrame.from_dict(raw_data)
    raw_df.head()
    raw_df.to_csv(output_dir + 'condensed_results.txt', sep='\t',index=False)


    for ds in datasets:
        analysis_id = '_'.join(['p','nomc','1','pearson','False',ds])
        df = raw_df[raw_df['analysis_id'] == analysis_id]

        # ensure white background per plot
        sns.set(font_scale=1.3)
        sns.set_style("ticks", {'font.family':'sans-serif','font.sans-serif':'Helvetica'})

        fig = plt.figure(figsize=(4,3))
        plt.title(ds_to_titles[ds], y=1.08)

        # generate subplot x ticks
        x_stats = ['p < 0.05', 'p > 0.05'] # latter are rho and tau

        # set labels
        labels = [['True Positive', 'reverse sign-True Positive',
                  'False Positive'], ['False Negative', 'True Negative']]
        labels = [['TP', 'rsTP', 'FP'], ['FN', 'TN']]


        # set colors # #228B22  #99ff99
        colors = [['#66b3ff','#ADD8E6','#ff9999'],['#228B22','#8064A2']]

        # iterate over for and rev dfs
        for d, name in enumerate(x_stats):

            # build bottom bar stack
            complete = np.zeros(1)
            for k, label in enumerate(labels[d]):
                # create bars
                plt.bar(d, df[label], bottom = complete, color=colors[d][k],
                        edgecolor='white', width=0.55, label=label, linewidth=0)#, tick_label=name)
                complete = np.add(complete, df[label])

        # dataset axes labels
        plt.xticks(ticks=[0,0.5,1], labels=['p < \u03B1','', 'p > \u03B1'])
        plt.ylabel('# of correlations')
        plt.tick_params(bottom=False)
        plt.tight_layout()
        sns.despine()

        fig.savefig(output_dir + 'barplots_dfcondensed_' + analysis_id + '.pdf')
        plt.close(fig)

        condensed_df = raw_df
        # p_fdr_1_rpearson_False_lungc

        if 'True' in cds:
            for p in params:
                for mc in mcs:
                    for fv in fvs:
                        for stat in statistics:
                            if stat[0] != 'r':
                                # create prefix analysis id, e.g. p_fdr_1
                                prefix = '_'.join([p, mc, fv, stat])

                                # ensure white background per plot with ticks
                                sns.set(font_scale=1.2)
                                sns.set_style("ticks", {'font.family':'sans-serif','font.sans-serif':'Helvetica'})

                                # create figure
                                fig, axarr = plt.subplots(nrows=1, ncols=3, figsize=(len(datasets)*3,4))

                                # Custom x axis
                                # plt.xlabel("Cook's D")

                                # generate subplot x ticks
                                x_stats = ['CUTIE','Cook\'s D']

                                # set labels
                                labels = ['True Positive', # 'reverse sign-True Positive',
                                          'False Positive', 'False Negative', 'True Negative']
                                labels = ['TP', 'FP', 'FN', 'TN']
                                # set colors # #228B22  #99ff99
                                colors = ['#66b3ff',# '#ADD8E6',
                                          '#ff9999','#228B22','#8064A2']

                                # iterate over datasets
                                for d, ds in enumerate(datasets):
                                    # get relevant subset of df
                                    # no then yes for Cooks D
                                    df = condensed_df[condensed_df['analysis_id'].isin([prefix + '_False_' + ds, prefix + '_True_' + ds])]


                                    # set number of bars (Cook D or not)
                                    r = range(2)

                                    # define subplot
                                    ax = plt.subplot(1, 3, d+1)
                                    plt.title(ds_to_titles[ds])
                                    # build bottom bar stack
                                    # fig = plt.figure(figsize=(8,4))
                                    complete = np.zeros(2)
                                    for k, label in enumerate(labels):
                                        # create bars
                                        plt.bar(r, df[label], bottom = complete, color=colors[k],
                                                edgecolor='white', width=0.85, label=label, linewidth=0)
                                        complete = np.add(complete, df[label])

                                    # subplot x ticks
                                    plt.xticks(r, x_stats)

                                    # remove axes
                                    ax.spines['right'].set_visible(False)
                                    ax.spines['top'].set_visible(False)

                                    # dataset x label
                                    # plt.xlabel('Cook\'s D')
                                    plt.tight_layout()

                                # Add a legend: not useful when rs-TP is negligible
                                # plt.legend(labels, loc='right', bbox_to_anchor=(1,1), ncol=1)

                                # ensures legend is not cropped
                                plt.tight_layout()

                                # save and close figure
                                fig.savefig(output_dir + 'Fig2CD_' + p + '_' + mc + '_' + fv + '_' + stat + '.pdf')
                                plt.close(fig)


if __name__ == "__main__":
    analyze_real()




