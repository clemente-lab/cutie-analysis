import pandas as pd
import numpy as np
import scipy.stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import glob
import os
import seaborn as sns
import click
import distutils.util

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.command(context_settings=CONTEXT_SETTINGS)
@click.version_option(version='0.1')

# Required arguments
@click.option('-fv', '--fold_value', type=str,
              help='fold value for criterion for p value change')
@click.option('-s', '--statistic', type=str,
              help='string denoting type of analysis')
@click.option('-p', '--param', type=str,
              help='string denoting type of param used')
@click.option('-c', '--corr_compare', type=str,
              help='boolean denoting whether performing cooksd or not')
@click.option('-cl', '--classes', type=str,
              help='types of input classes')
@click.option('-nse', '--n_seed', type=str,
              help='number of seeds used')
@click.option('-nsa', '--n_samp', type=str,
              help='number of samples used')
@click.option('-rn', '--rangestr', type=str,
              help='start stop and step of corr')
@click.option('-i', '--input_dir', type=click.Path(exists=True),
              help='input dir with .txt files of data')
@click.option('-o', '--output_dir', type=click.Path(exists=True),
              help='output dir to put config files')


def analyze_simulations(fold_value, statistic, param, corr_compare, classes,
    n_seed, n_samp, rangestr, input_dir, output_dir):
    '''
    Script for analysis of simulated data by CUTIE
    '''
    start, stop, step = [float(x) for x in rangestr.split(',')]

    # check if results df exists already
    if not os.path.exists(output_dir + 'sim_results_df.txt'):

        def parse_log(f, stat):
            lines = [l.strip() for l in f.readlines()]
            if stat != 'False':
                for l in lines:
                    if "number of correlations" in l:
                        n_corr = int(l.split(' ')[-1])
                    elif "initial_corr" in l:
                        initial_corr = int(l.split(' ')[-1])
                    elif "false correlations according to " + stat in l:
                        false_corr = int(l.split(' ')[-1])
                    elif "true correlations according to " + stat in l:
                        true_corr = int(l.split(' ')[-1])
                    elif "runtime" in l:
                        runtime = float(l.split(' ')[-1])
                rs_false = np.nan
                rs_true = np.nan

            else:
                # check if FDR correction defaulted
                for l in lines:
                    if "number of correlations" in l:
                        n_corr = int(l.split(' ')[-1])
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

            return n_corr, initial_corr, false_corr, true_corr, rs_false, rs_true, runtime


        df_dict = {}
        for p in param.split(','):
            df_dict[p] = {}
            for fv in fold_value.split(','):
                df_dict[p][fv] = {}
                for stat in statistic.split(','):
                    df_dict[p][fv][stat] = {}
                    for cc in corr_compare.split(','):
                        df_dict[p][fv][stat][cc] = {}
                        for seed in [str(x) for x in range(int(n_seed))]:
                            df_dict[p][fv][stat][cc][seed] = {}
                            for c in classes.split(','):
                                df_dict[p][fv][stat][cc][seed][c] = {}
                                for samp in n_samp.split(','):
                                    df_dict[p][fv][stat][cc][seed][c][samp] = {}
                                    for cor in ['{0:g}'.format(float(str(x))) for x in np.arange(start, stop+step, step)]:
                                        df_dict[p][fv][stat][cc][seed][c][samp][cor] = (np.nan, np.nan)


        file_dirs = glob.glob(input_dir + '*')
        missing = []
        done = 0
        failed = []

        # troubleshooting
        for f in file_dirs:
            subset_files = glob.glob(f + '/*.txt')
            subset_files.sort()
            try:
                # grab the most recent txt (log) file
                fn = subset_files[-1]
            except:
                print(f)

        for f in file_dirs:
            subset_files = glob.glob(f + '/*.txt')
            subset_files.sort()
            # grab the most recent txt (log) file
            fn = subset_files[-1]


            # get parts
            label = f.split('/')[-1]
            p, fv, stat, cc, seed, c, samp, cor = label.split('_')

            # when you grab the label split, if true, partition to the types in the args to the script
            if cc != 'False':
                for x in corr_compare.split(','):
                    if x != 'False':
                        #try:
                        with open(fn, 'r') as rf:
                            n_corr, initial_corr, false_corr, true_corr, rs_false, rs_true, runtime = parse_log(rf, stat=x)
                            df_dict[p][fv][stat][x][seed][c][samp][cor] = (true_corr, initial_corr)
                                #except:
                            #    failed.append(label)
                            #    print(label)

            # for CUTIE alone (AKA False)
            else:
                with open(fn, 'r') as rf:
                    n_corr, initial_corr, false_corr, true_corr, rs_false, rs_true, runtime = parse_log(rf, stat=cc)
                    df_dict[p][fv][stat][cc][seed][c][samp][cor] = (true_corr, initial_corr)

            done += 1
            if not subset_files:
                missing.append(f)

        missing.sort()
        # print([os.path.basename(x) for x in missing])
        analysis_ids = []
        ps = []
        fvs = []
        stats = []
        ccs = []
        seeds = []
        class_labs = []
        nsamps = []
        cors = []
        results = []
        for p in param.split(','):
            for fv in fold_value.split(','):
                for stat in statistic.split(','):
                    for cc in corr_compare.split(','):
                        for seed in [str(x) for x in range(int(n_seed))]:
                            for c in classes.split(','):
                                for samp in n_samp.split(','):
                                    for cor in ['{0:g}'.format(float(str(x))) for x in np.arange(start, stop+step, step)]:
                                        d = df_dict[p][fv][stat][cc][seed][c][samp][cor]
                                        # d = true corr, initial corr
                                        # if initial corr is 0, we don't add it to df
                                        if not np.isnan(d[0]):
                                            if d[1] == 1:
                                                analysis_ids.append('_'.join([p, fv, stat, cc, seed, c, samp, cor]))
                                                ps.append(p)
                                                fvs.append(fv)
                                                stats.append(stat)
                                                ccs.append(cc)
                                                seeds.append(seed)
                                                class_labs.append(c)
                                                nsamps.append(samp)
                                                cors.append(float(cor))
                                                results.append(float(d[0]))

        results_df = pd.DataFrame({'analysis_id': analysis_ids, 'parameter': ps,
                                   'fold_value': fvs, 'stat': stats, 'infln': ccs,
                                   'seed': seeds, 'class': class_labs,
                                   'sample_size': nsamps, 'corr_strength': cors,
                                   'indicator': results})

        results_df.to_csv(output_dir + 'sim_results_df.txt', sep='\t')
        # diagnostics
        print(len(missing),done,len(failed))
    else:
        # load in DF
        results_df = pd.read_csv(output_dir + 'sim_results_df.txt', sep='\t')
        # cast to string
        for col in ['fold_value', 'infln', 'sample_size']:
            results_df[col] = results_df[col].astype(str)
        # cast to float
        for col in ['corr_strength', 'indicator']:
            results_df[col] = results_df[col].astype(float)

    # corr_ticks = [0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, \
    # 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1]
    corr_ticks = ['0.0','','0.1','','0.2','','0.3','','0.4','','0.5','','0.6','',\
                  '0.7','','0.8','','0.9','','1.0']
    corr_ticks = ['{0:g}'.format(float(str(x))) for x in np.arange(start, stop+step, step)]
    for c in corr_ticks:
        if c[-1] != '0':
            c = ''

    infln_to_for = {
        'cookd': 'Cook\'s D (p < 0.05)',
        'dffits': 'DFFITS (p < 0.05)',
        'dsr': 'DSR (p < 0.05)'
    }

    infln_to_rev = {
        'cookd': 'Cook\'s D (p > 0.05)',
        'dffits': 'DFFITS (p > 0.05)',
        'dsr': 'DSR (p > 0.05)'

    }

    def new_label(row, stat):
        '''
        Relabels the statistic for the legend
        (1) If infln is True and the statistic is not pearson, don't make a plot for it
        (2) If infln is True and the statistic is pearson, then label with stat instead
        (3) If the statistic is TN/FN separation, label that line p > 0.05
        (4) Else the statistic is TP/FP separation, label that line p < 0.05
        '''
        # if an influence metric was used
        if row['infln'] != 'False':
            # ignore influence metrics for non pearson stats
            if row['stat'] in ['pearson','rpearson']:
                if row['stat'] == 'pearson':
                    return infln_to_for[stat]
                else:
                    return infln_to_rev[stat]
            else:
                return 'exclude'

        # follows this if Pearson only
        else:
            if row['stat'][0] == 'r':
                if stat != 'False':
                    return 'CUTIE (p > 0.05)'
                else:
                    return 'p > 0.05'
            else:
                if stat != 'False':
                    return 'CUTIE (p < 0.05)'
                else:
                    return 'p < 0.05'


    # grab statistics
    stat_pairs = []
    for v, w in zip(statistic.split(',')[::2], statistic.split(',')[1::2]):
        stat_pairs.append([v, w])

    # indiv plots
    for p in param.split(','):
        for fv in fold_value.split(','):
            for stat in stat_pairs:
                for cc in ['False']:
                    for c in classes.split(','):
                        for samp in n_samp.split(','):
                            # try statement exists because not all files will have an entry
                            try:
                                # subset dataframe
                                df = results_df[results_df['parameter'] == p]
                                df = df[df['fold_value'] == fv]
                                df = df[df['stat'].isin(stat)]
                                df = df[df['infln'] == cc]
                                df = df[df['class'] == c]
                                df = df[df['sample_size'] == samp]
                                df['Significance'] = df.apply(lambda row: new_label(row, cc),axis=1)
                                df = df.drop(['stat'], axis=1)

                                # for the FP simulation class, switch all 0s and 1s in pearson,
                                # but keep assignments in rpearson

                                # for the FN simulations, keep as is

                                # for CD simulations, keep as is ?
                                # if c == ['']
                                # df['indicator'] = 1 - df['indicator']

                                # set styles
                                sns.set(font_scale=1.4)
                                sns.set_style("ticks", {'font.family':'sans-serif','font.sans-serif':'Helvetica'})

                                # blue, red
                                colors = ['#4F81BD','#C0504D']
                                # stats = ['CUTIE (p < 0.05)', 'CUTIE (p > 0.05)']
                                stats = ['p < 0.05', 'p > 0.05']

                                if c == 'NP':
                                    ctitle = 'TP and TN'
                                else:
                                    ctitle = c

                                title = 'Power Curves for simulations of ' + ctitle + \
                                     '\n scatterplots using ' + stat[0].capitalize()  + ' and CUTIE'
                                plt.figure(figsize=(6,6))
                                ax = sns.pointplot(x="corr_strength", y="indicator", hue='Significance',data=df, ci=95,
                                    palette=sns.color_palette(colors), hue_order=stats)#, legend=False)
                                ax.set_title(title, fontsize=15)
                                plt.setp(ax.collections, alpha=.3) #for the markers
                                plt.setp(ax.lines, alpha=.3)
                                plt.ylim(-0.2, 1.2)

                                ax.set_ylabel('Proportion classified as True')
                                ax.set_xlabel('Correlation Strength')
                                ax.set_xticklabels(corr_ticks,rotation=0)
                                ax.set_yticklabels(['',0,0.2,0.4,0.6,0.8,1])

                                plt.tick_params(axis='both', which='both', top=False, right=False)
                                sns.despine()
                                plt.tight_layout()
                                df.to_csv(output_dir + '_'.join([p, fv, stat[0], cc, c, samp]) + '.csv', index=False)
                                plt.savefig(output_dir + '_'.join([p, fv, stat[0], cc, c, samp]) + '.pdf')
                                plt.close()
                            except:
                                 print('_'.join([p, fv, stat[0], cc, c, samp]))


    # cook D comparison
    for cc in corr_compare.split(','):
        for p in param.split(','):
            for fv in fold_value.split(','):
                for stat in [ ['pearson','rpearson'] ]:
                    for c in classes.split(','):
                        for samp in n_samp.split(','):
                            if cc != 'False':
                                #try:
                                # subset dataframe
                                df = results_df[results_df['parameter'] == p]
                                df = df[df['fold_value'] == fv]
                                df = df[df['stat'].isin(stat)]
                                df = df[df['class'] == c]
                                df = df[df['sample_size'] == samp]
                                df = df[df['infln'].isin(['False',cc])]
                                df['Method'] = df.apply(lambda row: new_label(row, cc),axis=1)
                                df = df[df['Method'] != 'exclude']
                                df = df.drop(['stat'], axis=1)

                                # generate plot
                                sns.set(font_scale=1.4)
                                sns.set_style("ticks", {'font.family':'sans-serif','font.sans-serif':'Helvetica'})


                                # green, blue, red
                                colors = ['#8064A2','#9BBB59','#4F81BD','#C0504D']
                                stats = [infln_to_rev[cc], infln_to_for[cc], 'CUTIE (p < 0.05)', 'CUTIE (p > 0.05)']

                                if c == 'NP':
                                    ctitle = 'TP and TN'
                                else:
                                    ctitle = c

                                title = 'Power Curves for simulations of ' + \
                                        c + '\n scatterplots using ' + stat[0].capitalize()

                                plt.figure(figsize=(6,6))
                                ax = sns.pointplot(x="corr_strength", y="indicator", hue='Method',data=df, ci=95,
                                    palette=sns.color_palette(colors), hue_order=stats)#, legend=False)
                                ax.set_title(title, fontsize=15)
                                plt.setp(ax.collections, alpha=.3) #for the markers
                                plt.setp(ax.lines, alpha=.3)
                                plt.ylim(-0.2,1.2)

                                ax.set_xticklabels(corr_ticks, rotation=0)
                                ax.set_yticklabels(['',0,0.2,0.4,0.6,0.8,1])
                                ax.set_ylabel('Proportion of Correlations classified as True')
                                ax.set_xlabel('Correlation Strength')

                                plt.tick_params(axis='both', which='both', top=False, right=False)
                                sns.despine()
                                plt.tight_layout()
                                df.to_csv(output_dir + '_'.join([p, fv, stat[0], cc, c, samp]) + '.csv', index=False)
                                plt.savefig(output_dir + '_'.join([p, fv, stat[0], cc, c, samp]) + '.pdf')
                                plt.close()
                                #except:
                                #    print(stat)
                                #   print(cc)


    print(results_df.head())

if __name__ == "__main__":
    analyze_simulations()




