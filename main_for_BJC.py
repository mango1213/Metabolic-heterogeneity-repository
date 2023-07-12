#%%
import os
from re import T
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None
from kaleido.scopes.plotly import PlotlyScope
import plotly.express as px
import plotly.graph_objects as go
from lifelines import CoxPHFitter
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

def load_data(folder_name_input):
    os.chdir(f'/Users/wangjun/Desktop/Script_wj/python/subpopulation_GC/{folder_name_input}')
    datas = pd.read_excel('1_df_clin_cluster_merged.xlsx')
    return datas

def preprocessing_GC_data(datas):
    if datas.columns.str.contains("Regions").any()==True:
        datas = datas.drop(["Regions"],axis=1)
    if datas.columns.str.contains("Number of spectra").any()==True:
        datas = datas.drop(["Number of spectra"],axis=1)
    if datas.columns.str.contains("Area [mm²]").any()==True:
        datas = datas.drop(["Area [mm²]"],axis=1)        
    datas = datas.set_index('Patient_ID')
    datas_piece = datas  
    datas_piece["my_calc_survtime_months"] = datas_piece.my_calc_survtime_months.astype(float)
    datas_piece["my_calc_censor"] = datas_piece.my_calc_censor.astype(float)
    return datas_piece

def Relabel_threshold (datas_piece,optimized_cutoff,set_threshold_b,cluster_id): 
    datas_plot_group = pd.DataFrame(columns=["survival","censor","group_id"])    
    for col_index, colum_name in enumerate(datas_piece.iloc[:,2:].columns.tolist()):
        criteria, label = colum_name.split("_")[-1].split(":")
        if set_threshold_b == True:
            datas_piece[colum_name] = datas_piece[colum_name].mask(datas_piece[colum_name] > float(optimized_cutoff)/100, int(label)) 
            datas_piece[colum_name] = datas_piece[colum_name].mask(datas_piece[colum_name] <= float(optimized_cutoff)/100, int(0)) 

        else:
            datas_piece[colum_name] = datas_piece[colum_name].mask(datas_piece[colum_name] > 1/float(criteria), int(label))
            datas_piece[colum_name] = datas_piece[colum_name].mask(datas_piece[colum_name] < 1/float(criteria), int(0))          
    if set_threshold_b==True:
        print(f'K={cluster_id}',f'threshold={optimized_cutoff}')
    else:
        criteria_cal=1/float(criteria)
        print(f'K={cluster_id}',f'threshold={criteria_cal}') 
    datas_piece_rm = datas_piece. loc[:, (datas_piece != 0).any(axis=0)] 
    datas_piece_rm = datas_piece_rm. loc[:,(datas_piece_rm != 0).sum() > 0] 

    for col_index, colum_name in enumerate(datas_piece_rm.iloc[:,2:].columns.tolist()):
        criteria, label = colum_name.split("_")[-1].split(":")   
        datas_pice_i = datas_piece_rm[datas_piece_rm[colum_name] == int(label)][['my_calc_survtime_months','my_calc_censor',colum_name]]        
        datas_pice_i.columns = datas_plot_group.columns
        datas_plot_group = pd.concat([datas_plot_group,datas_pice_i], axis=0, ignore_index=True)

    return datas_plot_group

def COX_AIC_within_cluster(datas_piece,threshold): 

    datas_piece = datas_piece.copy(deep=True)
    for col_index, colum_name in enumerate(datas_piece.iloc[:,2:].columns.tolist()):
        K, label = colum_name.split("_")[-1].split(":")
        datas_piece[colum_name] = datas_piece[colum_name].mask(datas_piece[colum_name] > float(threshold)/100, int(1)) 
        datas_piece[colum_name] = datas_piece[colum_name].mask(datas_piece[colum_name] < float(threshold)/100, int(0))
    datas_wo_null = datas_piece.loc[:, (datas_piece != 0).any(axis=0)] 
    datas_wo_null = datas_wo_null.loc[:,((datas_wo_null != 0).sum()> 4)]

    COX_re = CoxPHFitter(penalizer=0.001, l1_ratio=0.) 
    cph = COX_re.fit(datas_wo_null,'my_calc_survtime_months', event_col='my_calc_censor')
    AIC_within_threshold1=cph.AIC_partial_

    return AIC_within_threshold1

def AIC_heatmap(df_scaled,AIC_threshold,K):

    fig = px.imshow(df_scaled,
                    labels = dict(x="threshold [%]", y="K", color="Productivity")
                    )

    fig.add_trace(go.Scatter(
        x           = AIC_threshold,
        y           = K,
        legendgroup = "group2",
        name        = "thresholds based on the minimum AIC values",
        mode        = 'lines+markers',
        line        = dict(color= "White")
        ))
    fig.update_layout(
        title              = dict(text='AIC of Cox regresion models',x=0.5,y=0.98,xanchor='center',yanchor= 'middle'),
        font               = dict(family="Arial",size=32),
        yaxis              = dict(title='', titlefont_size=20,tickfont_size=32),
        xaxis              = dict(tickmode='linear',tick0 = 10, dtick = 6, titlefont_size=32,tickfont_size=32),
        coloraxis_colorbar = dict(title=''),
        showlegend         = True, 
        legend_orientation = 'h',
        legend             = dict(
                                x=.8, y=1.1,
                                font=dict(family='Arial',
                                size=20, color='black'),
                                bgcolor="LightSteelBlue"),
        plot_bgcolor       = 'white'
        )
    return fig

def AIC_heatmap_main(file,export1=False,export2=False):

    nb_of_clusters=['K=2','K=3','K=4','K=5','K=6','K=7','K=8','K=9','K=10']
    AIC_threshold = [16,10,26,4,8,26,20,24,24]
    x_axis_range = range(4,40,2)
    cluster_id_range = (np.arange(9)+2)     

    folder_name_output = f'{file}/AIC'  
    datas = load_data(folder_name_input)
    datas_piece= preprocessing_GC_data(datas)

    Dic_differ_K={}
    for cluster_id in cluster_id_range:  
        datas_piece_i = datas_piece.iloc[:,datas_piece.columns.str.contains(
            f"_{cluster_id}:|my_calc_survtime_months|my_calc_censor")]

        AIC_within_threshold = []; keys=[]; Dic={}
        for threshold in x_axis_range:
            AIC_within_threshold1 = COX_AIC_within_cluster(datas_piece_i,threshold)
            AIC_within_threshold.append(AIC_within_threshold1)
            keys.append(threshold)
        for m,k in zip(keys,AIC_within_threshold):
            Dic[m]=k
        Dic_differ_K['K={}'.format(cluster_id)]= Dic
    if export1 == True:
        Dic_differ_K = pd.DataFrame(Dic_differ_K)
        Dic_differ_K.to_excel('AIC_heatmap_input_6_51.xlsx')

    df_scaled = pd.DataFrame(Dic_differ_K)/200
    df_scaled_trans = df_scaled.transpose()
    fig             = AIC_heatmap(df_scaled_trans,AIC_threshold,nb_of_clusters)

    if export2 == True:
        os.chdir(f'/Users/wangjun/Desktop/Script_wj/python/subpopulation_GC/{folder_name_output}')
        if not os.path.exists("heatmap"):
            os.mkdir("heatmap")
        with open("./heatmap/AIC_heatmap.png", "wb") as f:   
           f.write(PlotlyScope().transform(fig, format="png",width=1600, height=1200))

def Relabel_threshold_for_ploting (datas_piece,threshold): 
    datas_plot_group = pd.DataFrame(columns=["survival","censor","group_id"])    
    for col_index, colum_name in enumerate(datas_piece.iloc[:,2:].columns.tolist()):
        criteria, label = colum_name.split("_")[-1].split(":")
        datas_piece[colum_name] = datas_piece[colum_name].mask(datas_piece[colum_name] > float(threshold)/100, int(label))
        datas_piece[colum_name] = datas_piece[colum_name].mask(datas_piece[colum_name] <= float(threshold)/100, int(0))
        datas_pice_i = datas_piece[datas_piece[colum_name] == int(label)][['my_calc_survtime_months','my_calc_censor',colum_name]]
        datas_pice_i.columns = datas_plot_group.columns
        datas_plot_group = pd.concat([datas_plot_group,datas_pice_i], axis=0, ignore_index=True)
    return datas_plot_group

def KM_subplot_2G (datas_piece_i, cluster_id, two_subgroups, cluster_id_index, nb_of_figure_x, nb_of_figure_y):
  
    ax = plt.subplot(nb_of_figure_x, nb_of_figure_y, cluster_id_index)   
    color = ["#dc143c", "#db7093",  "#0000ff","#868479", "#ff00ff", "#00ff7f", "#7b68ee", "#07ffff", '#da9d0c','#ff6347'] #color3 "#868479" wont be used here
    kmf = KaplanMeierFitter()
    for i in two_subgroups:     
        ix = (datas_piece_i['group_id']== i)
        
        kmf.fit(datas_piece_i['survival'][ix], event_observed = datas_piece_i['censor'][ix],
                        label = "subpopulation{}".format(i))
        bx = kmf.plot(title="K = {}".format(cluster_id), ax = ax, fontsize = 18,
                        ci_show = False, show_censors = True,              
                        color = color[i-1],
                        censor_styles = {'ms': 20, 'marker': '|',"mew": 4},
                        linewidth = 2.5)
    print(('K={}').format(cluster_id))
    print(('sub{}:n={}').format(two_subgroups[0],(datas_piece_i['group_id']== two_subgroups[0]).sum()))
    print(('sub{}:n={}').format(two_subgroups[1],(datas_piece_i['group_id']== two_subgroups[1]).sum()))

    bx.set_xlabel('Time(Mon)',fontsize=16)
    bx.set_ylabel('Survival Probability',fontsize=16)
    tl = ((bx.get_xlim()[1] - bx.get_xlim()[0])*0.010 + bx.get_xlim()[0],
          (bx.get_ylim()[1] - bx.get_ylim()[0])*0.15 + bx.get_ylim()[0])                
    bx.legend(fontsize=16,loc="upper right")
    bx.spines['right'].set_visible(False)
    bx.spines['top'].set_visible(False)  
    return bx  

def KM_subplot_all_comparison (datas_piece_i,cluster_id,cluster_id_index,nb_of_figure_x,nb_of_figure_y):

    ax = plt.subplot(nb_of_figure_x,nb_of_figure_y,cluster_id_index) 
    color = ["#dc143c", "#db7093", "#0000ff", "#868479", "#ff00ff", "#00ff7f", "#7b68ee", "#07ffff", '#da9d0c','#ff6347'] # 
    kmf = KaplanMeierFitter()
    for i in range(1,cluster_id+1):     
        ix = (datas_piece_i['group_id']== i)
        if datas_piece_i['survival'][ix].empty == False:
            kmf.fit(datas_piece_i['survival'][ix],event_observed=datas_piece_i['censor'][ix],label="subpopulation{}".format(i))
            bx = kmf.plot(  title="K = {}".format(cluster_id), 
                            ci_show       = False,                              show_censors = True,
                            color         = color[i-1],                         fontsize     = 16,
                            censor_styles = {'ms': 20, 'marker': '|',"mew": 4}, ax           = ax,
                            linewidth     = 2)

    bx.set_xlabel('Time(Mon)',fontsize=16)
    bx.set_ylabel('Survival Probability',fontsize=16)
    tl = ((bx.get_xlim()[1] - bx.get_xlim()[0])*0.010 + bx.get_xlim()[0],
          (bx.get_ylim()[1] - bx.get_ylim()[0])*0.15 + bx.get_ylim()[0])                
    bx.legend(fontsize=16,loc="upper right")
    bx.spines['right'].set_visible(False)
    bx.spines['top'].set_visible(False)  
    return bx

def KM_multiple_plots_pairwise_all_comparison_main(file, export = True):

    nb_of_figure_x = 2
    nb_of_figure_y = 4
    top_figure  = zip([9,9,9,9],[24,24,24,24],[[1,4],[1,5],[4,9],[5,9]],[1,2,3,4])
    down_figure = zip([9,9,9,9],[24,24,24,24],[5,6,7,8])

    folder_name_output = f'{file}/KM_plot'
    datas = load_data(folder_name_input)
    datas_piece= preprocessing_GC_data(datas)

    id = []
    # top 4 figures : two groups
    for cluster_id, threshold_top, two_subgroups,cluster_id_index in top_figure:
        datas_piece_i = datas_piece.iloc[:,datas_piece.columns.str.contains(
            f"_{cluster_id}:|my_calc_survtime_months|my_calc_censor")]   
        datas_piece_i = Relabel_threshold_for_ploting (datas_piece_i, threshold_top) 
        bx = KM_subplot_2G(datas_piece_i,cluster_id,two_subgroups,cluster_id_index,nb_of_figure_x,nb_of_figure_y) 

    # DOWN 4 figures: all groups
    for cluster_id, threshold_down, cluster_id_index in down_figure: 
        datas_piece_i = datas_piece.iloc[:,datas_piece.columns.str.contains(
            f"_{cluster_id}:|my_calc_survtime_months|my_calc_censor")] 
        datas_piece_i = Relabel_threshold_for_ploting (datas_piece_i, threshold_down)  
        bx = KM_subplot_all_comparison (datas_piece_i,cluster_id,cluster_id_index,nb_of_figure_x,nb_of_figure_y)
        id.append(cluster_id)
    if export == True:
        os.chdir(f'/Users/wangjun/Desktop/Script_wj/python/subpopulation_GC/{folder_name_output}')
        bx.get_figure().savefig(f"test_KM_K={id}.png")


if __name__ == '__main__':
    all = {'do_AIC_heatmap':  [],
           'do_KM_pairwise':  []}
    method = 'do_KM_pairwise'
    number_of_cluster = 9
    file = '6_3061_complete_mz'
    folder_name_input= f'{file}/input'
    if method == 'do_AIC_heatmap':
        AIC_heatmap_main(file)
    if method == 'do_KM_pairwise':
        KM_multiple_plots_pairwise_all_comparison_main(file)