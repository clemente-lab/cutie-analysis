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

        return defaulted, initial_corr, false_corr, true_corr, rs_false, rs_true, runtime


    headers = [
        'analysis_id',
        'parameter',
        'dataset',
        'statistic',
        'mc_used', #NEW
        'fold_value', # NEW
        'cooksd', #NEW
        'initial_corr',
        'true_corr(TP_FN)',
        'false_corr(FP_TN)',
        'rs_true_corr_TP_FN',
        'rs_false_corr_FP_TN',
        'runtime'
    ]

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
                                    defaulted, initial_corr, false_corr, \
                                        true_corr, rs_false, rs_true, runtime = parse_log(f,cd)

                                    new_row = pd.DataFrame([[analysis_id, p, d, s,
                                                            mc, fv, cd,
                                                            initial_corr, true_corr,
                                                            false_corr, rs_true,
                                                            rs_false, runtime]],
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

    # colnames = ['LungTranscriptomics', 'Micrometa', 'Microbiome', 'Gene Expression', 'WHO']
    colnames = ['LungCancer', 'LungTranscriptomics', 'Gene Expression', 'WHO']

    col_to_corr = {
        'LungTranscriptomics': 292 * 97, #depends on sum vs unsum
        'Micrometa': 83 * 897,
        'LungCancer': 748 * 747 / 2,
        'GeneExpression': 1000 * 999 / 2,
        'WHO': 354 * 353 / 2
    }

    # dists = ['lungtx', 'lungpt', 'lungc', 'hdac', 'who']
    dists = ['lungc','lungtx', 'hdac', 'who']

    dist_to_corr = {
        'lungtx': 292 * 97,
        'lungpt': 83 * 897,
        'lungc': 748 * 747 / 2,
        'hdac': 1000 * 999 / 2,
        'who': 354 * 353 / 2
    }
    stat_to_unicode = {
        'pearson': 'r',
        'spearman': '\u03C1', # rho
        'kendall': '\u03C4' # tau
    }

    results_df.to_csv(output_dir + 'real_results_df.txt', sep='\t')

    # populate indices and ids for the dfs
    for p in params:
        for fv in fvs:
            for mc in mcs:
                indices = []
                ids = []
                indices.append('_'.join(['pearson', 'cd', fv, mc, p]))
                indices.append('Pct initial corr')
                ids.append('_'.join([mc, fv, 'pearson', 'True', p]))
                for stat in stats:
                    indices.append('_'.join([stat, fv, mc]))
                    indices.append('Pct initial corr')
                    ids.append('_'.join([mc, fv, stat, 'False', p]))

                # populate  df
                df_array = []
                for i, (idstring, index) in enumerate(zip(ids, indices)):
                    row_fracs = []
                    mc, fv, s, cd, p = idstring.split('_')
                    for dist in dists:
                        row = results_df[(results_df['parameter'] == p) & (results_df['dataset'] == dist) & (results_df['statistic'] == s) \
                                     & (results_df['mc_used'] == mc) & (results_df['fold_value'] == fv) & (results_df['cooksd'] == cd)]
                        try:
                            row_fracs.append(float(row['true_corr(TP_FN)'] /row['initial_corr'].values)) # correctly id tp
                        except:
                            row_fracs.append(np.nan)
                            print('nan in row fracs')
                            print(dist, idstring)

                    df_array.append(row_fracs)

                    initial_sig_fracs = []
                    for dist in dists:
                        row = results_df[(results_df['parameter'] == p) & (results_df['dataset'] == dist) & (results_df['statistic'] == s) \
                                     & (results_df['mc_used'] == mc) & (results_df['fold_value'] == fv) & (results_df['cooksd'] == cd)]
                        # change number 249500 to n_corr depending on dataset
                        try:
                            initial_sig_fracs.append(float(row['initial_corr'] / dist_to_corr[dist]))
                        except:
                            initial_sig_fracs.append(np.nan)

                    df_array.append(initial_sig_fracs)

                pie_df = pd.DataFrame(data = df_array, index = indices, columns = colnames)
                pie_df = pie_df.rename_axis('Statistic')
                pie_df = pie_df.apply(pd.to_numeric).round(2)

                # parse the reverse sign shenanigans
                df_array = []

                # cut out the cookd parts
                rs_ids = ids[-len(stats):]
                rs_indices = indices[-2*len(stats):]
                for i, (idstring, index) in enumerate(zip(rs_ids, rs_indices)):
                    # stat = 'Pearson'
                    row_fracs = []
                    mc, fv, s, cd, p = idstring.split('_')
                    for dist in dists:
                        row = results_df[(results_df['parameter'] == p) & (results_df['dataset'] == dist) & (results_df['statistic'] == s) \
                                     & (results_df['mc_used'] == mc) & (results_df['fold_value'] == fv) & (results_df['cooksd'] == 'False')]
                        try:
                            row_fracs.append(float(row['rs_true_corr_TP_FN'] /row['initial_corr'].values)) # correctly id tp
                        except:
                            row_fracs.append(np.nan)
                            print('failed to parse rs')
                            print(dist, idstring)

                    df_array.append(row_fracs)

                    initial_sig_fracs = []
                    for dist in dists:
                        row = results_df[(results_df['parameter'] == p) & (results_df['dataset'] == dist) & (results_df['statistic'] == s) \
                                     & (results_df['mc_used'] == mc) & (results_df['fold_value'] == fv) & (results_df['cooksd'] == 'False')]
                        # change number 249500 to n_corr depending on dataset
                        try:
                            initial_sig_fracs.append(float(row['initial_corr'] / dist_to_corr[dist]))
                        except:
                            initial_sig_fracs.append(np.nan)

                    df_array.append(initial_sig_fracs)

                rs_df = pd.DataFrame(data = df_array, index = rs_indices, columns = colnames)
                rs_df = rs_df.rename_axis('Statistic')
                rs_df = rs_df.apply(pd.to_numeric).round(2)

                # currently the four dfs are
                # pie_df and rs_df
                # only pie_df has cookd info in it
                # the outer loop has mc and fv so when you save fig make sure to incl those

                # dictionary from which to get results for pie plots
                dd = {}

                # cut out micrometa dataset
                #pie_df = pie_df.drop(['Micrometa'],axis=1)
                nocd_pie_df = pie_df.iloc[2:,:]
                #rs_df = rs_df.drop(['Micrometa'],axis=1)
                sub_colnames = ['LungCancer', 'LungTranscriptomics', 'GeneExpression', 'WHO']

                # obtain indices without cook's D
                vals = list(nocd_pie_df.index.values)
                # skips by 2 (AKA every other)
                new_vals = vals[0::2]
                for v in new_vals:
                    dd[v] = {}

                for v in new_vals:
                    # v = 'pearson_1_fdr'
                    # check to make sure forward direction
                    if v.split('_')[0][0] != 'r':
                        dd[v]['rsTP'] = rs_df.loc[v,:].values
                    else:
                        dd[v]['rsFN'] = rs_df.loc[v,:].values


                for v in new_vals:
                    rows = nocd_pie_df.iloc[vals.index(v):vals.index(v)+2,:].values
                    if v.split('_')[0][0] != 'r':
                        dd[v]['TP'] = rows[0]
                        dd[v]['initial_sig'] = rows[1]
                    else:
                        dd[v]['FN'] = rows[0]
                        dd[v]['initial_insig'] = rows[1]


                for_vals = new_vals[::2]
                v_to_cd = {}
                # just get Cook's D
                # cd_val = list(pie_df.index.values)[0::2][0]

                # first two rows are cd
                rows = pie_df.iloc[0:2,:].values
                v_to_cd['TP'] = rows[0]
                v_to_cd['initial_sig'] = rows[1]

                # create figure
                # f, axarr = plt.subplots(len(for_vals) + 1,len(sub_colnames))
                # +1 to make room for the legend
                fig, axarr = plt.subplots(nrows=1, ncols=len(sub_colnames), # + 1,
                                          # sharex='col', sharey='row',
                                          figsize=(10,6))

                # Custom x axis
                plt.xlabel("Dataset")

                # generate subplot x ticks
                # extract statistic from first item in _ list
                # x = 'pearson_3_fdr'
                x_stats = [stat_to_unicode[x.split('_')[0]] for x in for_vals]

                # set labels
                labels = ['True Positive', 'reverse sign-True Positive',
                          'False Positive', 'False Negative', 'True Negative']

                # set colors
                colors = ['#66b3ff','#ADD8E6','#ff9999','#99ff99','#8064A2']

                # iterate over dataset
                first = True
                for d, name in enumerate(sub_colnames):

                    # iterate over statistic
                    stat_to_vals = defaultdict(list)

                    # hold values for stacked barplot
                    v_to_sizes = {}

                    for v, val in enumerate(for_vals):
                        # labels = ['TP', 'rsTP', 'FP', 'FN', 'rsFN', 'TN']

                        # TP is blue FP is red FN is green TN is purple
                        # for rs case
                        # reverse sign but still true FP is non reverse sign
                        TP = dd[val]['TP'][d]
                        rsTP = dd[val]['rsTP'][d]
                        P = dd[val]['initial_sig'][d]
                        FN = dd['r' + val]['FN'][d]
                        rsFN = dd['r' + val]['rsFN'][d]
                        N = dd['r' + val]['initial_insig'][d]
                        # sizes = [(TP - rsTP) * P, rsTP * P,(1-TP)*P, (FN - rsFN) * N, rsFN * N, (1-FN)*N]
                        sizes = [(TP - rsTP) * P, rsTP * P,(1-TP)*P, FN * N, (1-FN)*N]
                        v_to_sizes[val] = sizes

                    # create df
                    raw_data = defaultdict(list)
                    for j, l in enumerate(labels):
                        for v in for_vals:
                            raw_data[l].append(v_to_sizes[v][j])


                    df = pd.DataFrame(raw_data)

                    # set number of bars (# of statistics)
                    r = range(len(for_vals))

                    # define subplot
                    # plt.subplot(1, len(colnames), d+1)
                    ax = plt.subplot(1, len(colnames), d+1)

                    # ensure white background per plot
                    sns.set_style('white')

                    # build bottom bar stack
                    # fig = plt.figure(figsize=(8,4))
                    complete = np.zeros(len(for_vals))
                    for k, label in enumerate(labels):
                        # create bars
                        plt.bar(r, df[label], bottom = complete, color=colors[k],
                                edgecolor='white', width=0.85, label=label)
                        complete = np.add(complete, df[label])

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


                # Add a legend
                plt.legend(labels, loc='right', bbox_to_anchor=(1,1), ncol=1)

                # ensures legend is not cropped
                plt.tight_layout()

                # save and close figure
                fig.savefig(output_dir + 'barplots_dfreal_combined_' + p + '_' + mc + '_' + fv + '.pdf')
                plt.close(fig)


if __name__ == "__main__":
    analyze_simulations_real()




