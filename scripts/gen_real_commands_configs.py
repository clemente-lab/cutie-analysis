import glob
import os
import click
import numpy as np
from cutie import parse

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
              help='string denoting parameter used')
@click.option('-d', '--datasets', type=str,
              help='string denoting datasets used')
@click.option('-c', '--corr_compare', type=str,
              help='boolean denoting whether performing cooksd or not')
@click.option('-cf', '--cutie_fp', type=click.Path(exists=True),
              help='path of cutie script executable')
@click.option('-w', '--working_dir', type=click.Path(exists=True),
              help='working dir to save results')
@click.option('-o', '--output_dir', type=click.Path(exists=True),
              help='output dir to put config files')

def gen_commands_configs(fold_value, statistic, multi_corr, param, datasets,
                         corr_compare, cutie_fp, working_dir, output_dir):
    data_to_params = {
        'hdac': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/data/HDAC_data/GSE15222_series_matrix_x100_del62.txt',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/data/HDAC_data/GSE15222_series_matrix_x100_del62.txt',
            'delimiter1': '\\t',
            'delimiter2': '\\t',
            'f1type': 'untidy',
            'f2type': 'untidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'True',
            'alpha': '0.05'},
        'lungtx': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/data/lungtx_data/otu_table_L6_filt1e3.txt',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/data/lungtx_data/Genes.KEGG.L3.add_counts.txt',
            'delimiter1': '\\t',
            'delimiter2': '\\t',
            'f1type': 'untidy',
            'f2type': 'untidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'False',
            'alpha': '0.05'},
        'lungc': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/data/pre_sparcc_MSQ/otu_table.MSQ34_L6.txt',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/data/pre_sparcc_MSQ/otu_table.MSQ34_L6.txt',
            'delimiter1': '\\t',
            'delimiter2': '\\t',
            'f1type': 'untidy',
            'f2type': 'untidy',
            'skip1': '1',
            'skip2': '1',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'True',
            'alpha': '0.05'},
        'lungpt': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/data/lungpt_data/otu_table_MultiO_merged___L6.txt',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/data/lungpt_data/Mapping.Pneumotype.Multiomics.RL.NYU.w_metabolites.w_inflamm.txt',
            'delimiter1': '\\t',
            'delimiter2': '\\t',
            'f1type': 'untidy',
            'f2type': 'tidy',
            'skip1': '1',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '16',
            'endcol2': '99',
            'paired': 'False',
            'alpha': '0.05'},
        'who': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/data/MINE_data/WHOfix.txt',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/data/MINE_data/WHOfix.txt',
            'delimiter1': '\\t',
            'delimiter2': '\\t',
            'f1type': 'tidy',
            'f2type': 'tidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '2',
            'endcol1': '356',
            'startcol2': '2',
            'endcol2': '356',
            'paired': 'True',
            'alpha': '0.05'},
        'baseball': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/MINE/Baseballfix.txt',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/MINE/Baseballfix.txt',
            'delimiter1': '\\t',
            'delimiter2': '\\t',
            'f1type': 'tidy',
            'f2type': 'tidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'True',
            'alpha': '0.05'},
        'ad0': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/mennonites/atopic_dz0.csv',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/mennonites/atopic_dz0.csv',
            'delimiter1': ',',
            'delimiter2': ',',
            'f1type': 'tidy',
            'f2type': 'tidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'True',
            'alpha': '0.05'},
        'ad1': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/mennonites/atopic_dz1.csv',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/mennonites/atopic_dz1.csv',
            'delimiter1': ',',
            'delimiter2': ',',
            'f1type': 'tidy',
            'f2type': 'tidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'True',
            'alpha': '0.05'},
        'oom': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/mennonites/OOM.csv',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/mennonites/OOM.csv',
            'delimiter1': ',',
            'delimiter2': ',',
            'f1type': 'tidy',
            'f2type': 'tidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'True',
            'alpha': '0.05'},
        'roc': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/mennonites/ROC.csv',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/mennonites/ROC.csv',
            'delimiter1': ',',
            'delimiter2': ',',
            'f1type': 'tidy',
            'f2type': 'tidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'True',
            'alpha': '0.05'},
        'mennonites': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/mennonites/df_mennonites.txt',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/mennonites/df_mennonites.txt',
            'delimiter1': '\\t',
            'delimiter2': '\\t',
            'f1type': 'tidy',
            'f2type': 'tidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'True',
            'alpha': '0.05'},
        'covid': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/df_covid.txt',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/df_covid.txt',
            'delimiter1': '\\t',
            'delimiter2': '\\t',
            'f1type': 'tidy',
            'f2type': 'tidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'True',
            'alpha': '0.05'},
        'airplane': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/data/MINE_data/2008_data.txt',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/data/MINE_data/2008_data.txt',
            'delimiter1': '\\t',
            'delimiter2': '\\t',
            'f1type': 'tidy',
            'f2type': 'tidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'True',
            'alpha': '0.05'},
        'ici': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/df_pre.txt',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/df_ici.txt',
            'delimiter1': '\\t',
            'delimiter2': '\\t',
            'f1type': 'tidy',
            'f2type': 'tidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'False',
            'alpha': '0.05'},
        'liverf': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/liver/df_liver_female_500.txt',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/liver/df_liver_female_500.txt',
            'delimiter1': '\\t',
            'delimiter2': '\\t',
            'f1type': 'untidy',
            'f2type': 'untidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'True',
            'alpha': '0.05'},
        'liverm': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/liver/df_liver_male_500.txt',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/liver/df_liver_male_500.txt',
            'delimiter1': '\\t',
            'delimiter2': '\\t',
            'f1type': 'untidy',
            'f2type': 'untidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'True',
            'alpha': '0.05'},
        'micro': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/MINE/Microbiome_fix_500.csv',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/MINE/Microbiome_fix_500.csv',
            'delimiter1': ',',
            'delimiter2': ',',
            'f1type': 'tidy',
            'f2type': 'tidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'True',
            'alpha': '0.05'},
        'hgoral': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/df_oral.txt',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/df_oral.txt',
            'delimiter1': '\\t',
            'delimiter2': '\\t',
            'f1type': 'tidy',
            'f2type': 'tidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'True',
            'alpha': '0.05'},
        'crc': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/df_CRC_otu.txt',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/df_CRC_cyto.txt',
            'delimiter1': '\\t',
            'delimiter2': '\\t',
            'f1type': 'tidy',
            'f2type': 'tidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'False',
            'alpha': '0.05'},
        'ibd': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/df_nat_otu.txt',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/df_nat_meta.txt',
            'delimiter1': '\\t',
            'delimiter2': '\\t',
            'f1type': 'tidy',
            'f2type': 'tidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'False',
            'alpha': '0.0047'},
        'cell': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/df_cell_Al.txt',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/df_cell_Al.txt',
            'delimiter1': '\\t',
            'delimiter2': '\\t',
            'f1type': 'tidy',
            'f2type': 'tidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'True',
            'alpha': '0.05'},
        'nc': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/natcom_fix.txt',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/natcom_fix.txt',
            'delimiter1': '\\t',
            'delimiter2': '\\t',
            'f1type': 'tidy',
            'f2type': 'tidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'True',
            'alpha': '0.05'},
        'plos': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/plos_fungi.txt',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/plos_bact.txt',
            'delimiter1': '\\t',
            'delimiter2': '\\t',
            'f1type': 'tidy',
            'f2type': 'tidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'False',
            'alpha': '0.05'},
        'ca': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/ca_otu.csv',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/ca_plasma.csv',
            'delimiter1': ',',
            'delimiter2': ',',
            'f1type': 'tidy',
            'f2type': 'tidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'False',
            'alpha': '0.05'},
        'statin': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/bmis/df_combined.csv',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/bmis/df_combined.csv',
            'delimiter1': ',',
            'delimiter2': ',',
            'f1type': 'tidy',
            'f2type': 'tidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'True',
            'alpha': '0.01128'},
        'spatial': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/spatial/df_f1.csv',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/spatial/df_f1.csv',
            'delimiter1': ',',
            'delimiter2': ',',
            'f1type': 'tidy',
            'f2type': 'tidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'True',
            'alpha': '0.05'},
        'livermfull': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/liver/df_liver_male.txt',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/liver/df_liver_male.txt',
            'delimiter1': '\\t',
            'delimiter2': '\\t',
            'f1type': 'untidy',
            'f2type': 'untidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'False',
            'alpha': '0.05',
            'njobs': 1000},
        'liverffull': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/liver/df_liver_female.txt',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/liver/df_liver_female.txt',
            'delimiter1': '\\t',
            'delimiter2': '\\t',
            'f1type': 'untidy',
            'f2type': 'untidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'False',
            'alpha': '0.05',
            'njobs': 1000},
        'microfull': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/MINE/Microbiome_fix.csv',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/cutie_exp/inputs/MINE/Microbiome_fix.csv',
            'delimiter1': ',',
            'delimiter2': ',',
            'f1type': 'tidy',
            'f2type': 'tidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'True',
            'alpha': '0.05',
            'njobs': 1000},
        'hdacfull': {
            'samp_var1_fp': '/sc/arion/projects/clemej05a/kevin/data/HDAC_data/GSE15222_series_matrix_full_del62.txt',
            'samp_var2_fp': '/sc/arion/projects/clemej05a/kevin/data/HDAC_data/GSE15222_series_matrix_full_del62.txt',
            'delimiter1': '\\t',
            'delimiter2': '\\t',
            'f1type': 'untidy',
            'f2type': 'untidy',
            'skip1': '0',
            'skip2': '0',
            'startcol1': '-1',
            'endcol1': '-1',
            'startcol2': '-1',
            'endcol2': '-1',
            'paired': 'True',
            'alpha': '0.05',
            'njobs': 1000},

        }
        # liverffull, microfull, hdacfull
        # endcol startcol check: lungpt and WHO

    fv = fold_value
    # files = glob.glob(input_dir + '*.txt')
    datasets = datasets.split(',')
    # datasets = ['hdac','lungc','lungpt','who','tx']
    for data in datasets:
        param_to_str = data_to_params[data]

        # example fid: p_nomc_1_mine_False_lungtx
        f_id = '_'.join([param, multi_corr, fv, statistic, corr_compare, data])

        ftype, samp_var_fp, startcol, endcol, delimiter, skip = param_to_str['f1type'], \
            param_to_str['samp_var1_fp'], int(param_to_str['startcol1']), \
            int(param_to_str['endcol1']), param_to_str['delimiter1'], int(param_to_str['skip1'])
        samp_ids, var_names, samp_var_df, n_var, n_samp = parse.parse_input(
            ftype, samp_var_fp, startcol, endcol, delimiter, skip)

        try:
            njobs = param_to_str['njobs']
        except:
            njobs = 1

        # create column tuples
        if njobs > 1:
            # create subtypes
            # samp_var_df is always in tidy format and has already been iloc'd
            dfs = np.array_split(samp_var_df, njobs, axis=1)
            vals = [df.shape[1] for df in dfs]

            col_tuples = [(0,0)]

            indices = [0]
            indices.extend(vals)
            for i in range(len(indices)-1):
                t = []
                prev = col_tuples[i]
                t.append(prev[1])
                t.append(indices[i+1] + prev[1])
                col_tuples.append(t)

            # get rid of 0,0  placeholder
            col_tuples.pop(0)
        else:
            col_tuples = [[param_to_str['endcol1'],param_to_str['endcol2']]]

        for i in range(njobs):
            # sub fid
            if njobs > 1:
                fid = f_id + '_' + str(i)
            else:
                fid = f_id

            # output_dir = '/sc/arion/projects/clemej05a/kevin/data/real_data_analysis/'
            # out_dir = output_dir + f_id + '/'
            out_dir = output_dir + fid + '/'
            try:
                os.makedirs(out_dir)
            except:
                pass
            # working_dir = '/sc/hydra/scratch/buk02/real_data_analysis/'
            working_outdir = working_dir + fid + '/'
            try:
                os.makedirs(working_outdir)
            except:
                pass

            with open(out_dir + 'config_' + fid + '.txt','w') as f:
                f.write('[input]')
                f.write('\n')
                f.write('samp_var1_fp: ' + param_to_str['samp_var1_fp'])
                f.write('\n')
                f.write('delimiter1: ' + param_to_str['delimiter2'])
                f.write('\n')
                f.write('samp_var2_fp: ' + param_to_str['samp_var2_fp'])
                f.write('\n')
                f.write('delimiter2: ' + param_to_str['delimiter2'])
                f.write('\n')
                f.write('f1type: ' + param_to_str['f1type'])
                f.write('\n')
                f.write('f2type: ' + param_to_str['f2type'])
                f.write('\n')
                f.write('skip1: ' + param_to_str['skip1'])
                f.write('\n')
                f.write('skip2: ' + param_to_str['skip2'])
                f.write('\n')
                f.write('startcol1: ' + param_to_str['startcol1'])
                f.write('\n')
                f.write('endcol1: ' + param_to_str['endcol1'])
                f.write('\n')
                f.write('startcol2: ' + str(col_tuples[i][0]))
                f.write('\n')
                f.write('endcol2: ' + str(col_tuples[i][1]))
                f.write('\n')
                f.write('paired: ' + param_to_str['paired'])
                f.write('\n')
                f.write('\n')
                f.write('[output]')
                f.write('\n')
                f.write('working_dir: ' + working_outdir)
                f.write('\n')
                f.write('overwrite: True')
                f.write('\n')
                f.write('\n')
                f.write('[stats]')
                f.write('\n')
                f.write('param: ' + param)
                f.write('\n')
                f.write('statistic: ' + statistic)
                f.write('\n')
                f.write('resample_k: 1')
                f.write('\n')
                f.write('alpha: ' + param_to_str['alpha'])
                f.write('\n')
                f.write('mc: ' + multi_corr)
                f.write('\n')
                f.write('fold: True')
                f.write('\n')
                f.write('fold_value: ' + fv)
                f.write('\n')
                f.write('corr_compare: ' + corr_compare)
                f.write('\n')
                f.write('\n')
                f.write('[graph]')
                f.write('\n')
                f.write('graph_bound: 30')
                f.write('\n')
                f.write('fix_axis: False')

            with open(out_dir + 'commands_' + fid + '.txt','w') as f:
                f.write('export PYTHONPATH=$PYTHONPATH:/hpc/users/buk02/tools/sandbox/lib/python3.7/site-packages/ && python ' + \
                        cutie_fp + ' -i ' + out_dir + 'config_' + fid + '.txt')

if __name__ == "__main__":
    gen_commands_configs()
