import pandas as pd
import numpy as np
import plotly.express as px
import matplotlib.pyplot as plt

sorted_array = [
    # Clustered del
    'clustered_del_1-10Kb_very early',
    'clustered_del_1-10Kb_early',
    'clustered_del_1-10Kb_medium',
    'clustered_del_1-10Kb_late',
    'clustered_del_1-10Kb_very late',

    'clustered_del_10-100Kb_very early',
    'clustered_del_10-100Kb_early',
    'clustered_del_10-100Kb_medium',
    'clustered_del_10-100Kb_late',
    'clustered_del_10-100Kb_very late',

    'clustered_del_100Kb-1Mb_very early',
    'clustered_del_100Kb-1Mb_early',
    'clustered_del_100Kb-1Mb_medium',
    'clustered_del_100Kb-1Mb_late',
    'clustered_del_100Kb-1Mb_very late',

    'clustered_del_1Mb-10Mb_very early',
    'clustered_del_1Mb-10Mb_early',
    'clustered_del_1Mb-10Mb_medium',
    'clustered_del_1Mb-10Mb_late',
    'clustered_del_1Mb-10Mb_very late',

    'clustered_del_>10Mb_very early',
    'clustered_del_>10Mb_early',
    'clustered_del_>10Mb_medium',
    'clustered_del_>10Mb_late',
    'clustered_del_>10Mb_very late',

    # Clustered tds
    'clustered_tds_1-10Kb_very early',
    'clustered_tds_1-10Kb_early',
    'clustered_tds_1-10Kb_medium',
    'clustered_tds_1-10Kb_late',
    'clustered_tds_1-10Kb_very late',

    'clustered_tds_10-100Kb_very early',
    'clustered_tds_10-100Kb_early',
    'clustered_tds_10-100Kb_medium',
    'clustered_tds_10-100Kb_late',
    'clustered_tds_10-100Kb_very late',

    'clustered_tds_100Kb-1Mb_very early',
    'clustered_tds_100Kb-1Mb_early',
    'clustered_tds_100Kb-1Mb_medium',
    'clustered_tds_100Kb-1Mb_late',
    'clustered_tds_100Kb-1Mb_very late',

    'clustered_tds_1Mb-10Mb_very early',
    'clustered_tds_1Mb-10Mb_early',
    'clustered_tds_1Mb-10Mb_medium',
    'clustered_tds_1Mb-10Mb_late',
    'clustered_tds_1Mb-10Mb_very late',

    'clustered_tds_>10Mb_very early',
    'clustered_tds_>10Mb_early',
    'clustered_tds_>10Mb_medium',
    'clustered_tds_>10Mb_late',
    'clustered_tds_>10Mb_very late',

    # Clustered inv
    'clustered_inv_1-10Kb_very early',
    'clustered_inv_1-10Kb_early',
    'clustered_inv_1-10Kb_medium',
    'clustered_inv_1-10Kb_late',
    'clustered_inv_1-10Kb_very late',

    'clustered_inv_10-100Kb_very early',
    'clustered_inv_10-100Kb_early',
    'clustered_inv_10-100Kb_medium',
    'clustered_inv_10-100Kb_late',
    'clustered_inv_10-100Kb_very late',

    'clustered_inv_100Kb-1Mb_very early',
    'clustered_inv_100Kb-1Mb_early',
    'clustered_inv_100Kb-1Mb_medium',
    'clustered_inv_100Kb-1Mb_late',
    'clustered_inv_100Kb-1Mb_very late',

    'clustered_inv_1Mb-10Mb_very early',
    'clustered_inv_1Mb-10Mb_early',
    'clustered_inv_1Mb-10Mb_medium',
    'clustered_inv_1Mb-10Mb_late',
    'clustered_inv_1Mb-10Mb_very late',

    'clustered_inv_>10Mb_very early',
    'clustered_inv_>10Mb_early',
    'clustered_inv_>10Mb_medium',
    'clustered_inv_>10Mb_late',
    'clustered_inv_>10Mb_very late',

    # Clustered trans
    'clustered_trans_very early',
    'clustered_trans_early',
    'clustered_trans_medium',
    'clustered_trans_late',
    'clustered_trans_very late',

    # Non-clustered del
    'non-clustered_del_1-10Kb_very early',
    'non-clustered_del_1-10Kb_early',
    'non-clustered_del_1-10Kb_medium',
    'non-clustered_del_1-10Kb_late',
    'non-clustered_del_1-10Kb_very late',

    'non-clustered_del_10-100Kb_very early',
    'non-clustered_del_10-100Kb_early',
    'non-clustered_del_10-100Kb_medium',
    'non-clustered_del_10-100Kb_late',
    'non-clustered_del_10-100Kb_very late',

    'non-clustered_del_100Kb-1Mb_very early',
    'non-clustered_del_100Kb-1Mb_early',
    'non-clustered_del_100Kb-1Mb_medium',
    'non-clustered_del_100Kb-1Mb_late',
    'non-clustered_del_100Kb-1Mb_very late',

    'non-clustered_del_1Mb-10Mb_very early',
    'non-clustered_del_1Mb-10Mb_early',
    'non-clustered_del_1Mb-10Mb_medium',
    'non-clustered_del_1Mb-10Mb_late',
    'non-clustered_del_1Mb-10Mb_very late',

    'non-clustered_del_>10Mb_very early',
    'non-clustered_del_>10Mb_early',
    'non-clustered_del_>10Mb_medium',
    'non-clustered_del_>10Mb_late',
    'non-clustered_del_>10Mb_very late',

    # Clustered tds
    'non-clustered_tds_1-10Kb_very early',
    'non-clustered_tds_1-10Kb_early',
    'non-clustered_tds_1-10Kb_medium',
    'non-clustered_tds_1-10Kb_late',
    'non-clustered_tds_1-10Kb_very late',

    'non-clustered_tds_10-100Kb_very early',
    'non-clustered_tds_10-100Kb_early',
    'non-clustered_tds_10-100Kb_medium',
    'non-clustered_tds_10-100Kb_late',
    'non-clustered_tds_10-100Kb_very late',

    'non-clustered_tds_100Kb-1Mb_very early',
    'non-clustered_tds_100Kb-1Mb_early',
    'non-clustered_tds_100Kb-1Mb_medium',
    'non-clustered_tds_100Kb-1Mb_late',
    'non-clustered_tds_100Kb-1Mb_very late',

    'non-clustered_tds_1Mb-10Mb_very early',
    'non-clustered_tds_1Mb-10Mb_early',
    'non-clustered_tds_1Mb-10Mb_medium',
    'non-clustered_tds_1Mb-10Mb_late',
    'non-clustered_tds_1Mb-10Mb_very late',

    'non-clustered_tds_>10Mb_very early',
    'non-clustered_tds_>10Mb_early',
    'non-clustered_tds_>10Mb_medium',
    'non-clustered_tds_>10Mb_late',
    'non-clustered_tds_>10Mb_very late',

    # Clustered inv
    'non-clustered_inv_1-10Kb_very early',
    'non-clustered_inv_1-10Kb_early',
    'non-clustered_inv_1-10Kb_medium',
    'non-clustered_inv_1-10Kb_late',
    'non-clustered_inv_1-10Kb_very late',

    'non-clustered_inv_10-100Kb_very early',
    'non-clustered_inv_10-100Kb_early',
    'non-clustered_inv_10-100Kb_medium',
    'non-clustered_inv_10-100Kb_late',
    'non-clustered_inv_10-100Kb_very late',

    'non-clustered_inv_100Kb-1Mb_very early',
    'non-clustered_inv_100Kb-1Mb_early',
    'non-clustered_inv_100Kb-1Mb_medium',
    'non-clustered_inv_100Kb-1Mb_late',
    'non-clustered_inv_100Kb-1Mb_very late',

    'non-clustered_inv_1Mb-10Mb_very early',
    'non-clustered_inv_1Mb-10Mb_early',
    'non-clustered_inv_1Mb-10Mb_medium',
    'non-clustered_inv_1Mb-10Mb_late',
    'non-clustered_inv_1Mb-10Mb_very late',

    'non-clustered_inv_>10Mb_very early',
    'non-clustered_inv_>10Mb_early',
    'non-clustered_inv_>10Mb_medium',
    'non-clustered_inv_>10Mb_late',
    'non-clustered_inv_>10Mb_very late',

    # Non-clustered trans
    'non-clustered_trans_very early',
    'non-clustered_trans_early',
    'non-clustered_trans_medium',
    'non-clustered_trans_late',
    'non-clustered_trans_very late',
]

color_map = {
    # Deletions - shades of red
    'Clustered Deletion (very early)': '#c61033',  # Light red
    'Clustered Deletion (early)': '#c61033',  # Medium red
    'Clustered Deletion (medium)': '#c61033',  # Dark red
    'Clustered Deletion (late)': '#c61033',  # Bright red
    'Clustered Deletion (very late)': '#c61033',  # Deep red
    'Non-clustered Deletion (very early)': '#c61033',  # Light red
    'Non-clustered Deletion (early)': '#c61033',  # Medium red
    'Non-clustered Deletion (medium)': '#c61033',  # Dark red
    'Non-clustered Deletion (late)': '#c61033',  # Bright red
    'Non-clustered Deletion (very late)': '#c61033',  # Deep red
    
    # Inversions - shades of blue
    'Clustered Inversion (very early)': '#3b74ac',  # Light blue
    'Clustered Inversion (early)': '#3b74ac',  # Medium blue
    'Clustered Inversion (medium)': '#3b74ac',  # Dark blue
    'Clustered Inversion (late)': '#3b74ac',  # Bright blue
    'Clustered Inversion (very late)': '#3b74ac',  # Deep blue
    'Non-clustered Inversion (very early)': '#3b74ac',  # Light blue
    'Non-clustered Inversion (early)': '#3b74ac',  # Medium blue
    'Non-clustered Inversion (medium)': '#3b74ac',  # Dark blue
    'Non-clustered Inversion (late)': '#3b74ac',  # Bright blue
    'Non-clustered Inversion (very late)': '#3b74ac',  # Deep blue
    
    # Tandem Duplications - shades of green
    'Clustered Tandem Duplication (very early)': '#44a23c',  # Light green
    'Clustered Tandem Duplication (early)': '#44a23c',  # Medium green
    'Clustered Tandem Duplication (medium)': '#44a23c',  # Dark green
    'Clustered Tandem Duplication (late)': '#44a23c',  # Bright green
    'Clustered Tandem Duplication (very late)': '#44a23c',  # Deep green
    'Non-clustered Tandem Duplication (very early)': '#44a23c',  # Light green
    'Non-clustered Tandem Duplication (early)': '#44a23c',  # Medium green
    'Non-clustered Tandem Duplication (medium)': '#44a23c',  # Dark green
    'Non-clustered Tandem Duplication (late)': '#44a23c',  # Bright green
    'Non-clustered Tandem Duplication (very late)': '#44a23c',  # Deep green
    
    # Translocations - shades of purple
    'Clustered Translocation (very early)': '#833792',  # Light purple
    'Clustered Translocation (early)': '#833792',  # Medium purple
    'Clustered Translocation (medium)': '#833792',  # Dark purple
    'Clustered Translocation (late)': '#833792',  # Bright purple
    'Clustered Translocation (very late)': '#833792',  # Deep purple
    'Non-clustered Translocation (very early)': '#833792',  # Light purple
    'Non-clustered Translocation (early)': '#833792',  # Medium purple
    'Non-clustered Translocation (medium)': '#833792',  # Dark purple
    'Non-clustered Translocation (late)': '#833792',  # Bright purple
    'Non-clustered Translocation (very late)': '#833792',  # Deep purple
}

colDict = {}
colDict['del'] = '#c61033'
colDict['tds'] = '#44a23c'
colDict['inv'] = '#3b74ac'
colDict['trans'] = '#833792'

def plot_mut_spectrum(mutation_vector, sorted_array, color_map, title, output_dir, show=True):
    def categorize_mutation(mutation):
        cluster = "Clustered" if mutation.startswith("clustered_") else "Non-clustered"
        if "_del_" in mutation:
            mtype = "Deletion"
        elif "_inv_" in mutation:
            mtype = "Inversion"
        elif "_tds_" in mutation:
            mtype = "Tandem Duplication"
        elif "_trans_" in mutation:
            mtype = "Translocation"
        else:
            mtype = "Other"
        subcat = mutation.split("_")[2] if mtype == "Translocation" else mutation.split("_")[3]
        return f"{cluster} {mtype} ({subcat})"

    df = mutation_vector.copy()
    df["Group"] = df["MutationType"].apply(categorize_mutation)
    df["Color"] = df["Group"].map(color_map)
    df["MutationType"] = pd.Categorical(df["MutationType"], categories=sorted_array, ordered=True)
    df.sort_values("MutationType", inplace=True)
    df.reset_index(drop=True, inplace=True)

    fig = px.bar(
        df, x="MutationType", y="Total Count", color="Group", title=title,
        labels={"Total Count": "Mutation Count", "MutationType": "Mutation Type"},
        color_discrete_map=color_map, text_auto=False,
        category_orders={"MutationType": sorted_array}
    )
    fig.update_traces(marker_line_width=0)
    fig.update_layout(title_x=0.5)

    y_max = df["Total Count"].max()
    for i in range(0, 161, 5):
        fig.add_vline(
            x=i - 0.5, line_dash="dot", line_color="grey", opacity=1, line_width=0.3
        )

    extra_bars = [
        (0, 25, "Deletions", '#c61033'),
        (25, 50, "Duplications", '#44a23c'),
        (50, 75, "Inversions", '#3b74ac'),
        (75, 80, "Trans", '#833792'),
        (80, 105, "Deletions", '#c61033'),
        (105, 130, "Duplications", '#44a23c'),
        (130, 155, "Inversions", '#3b74ac'),
        (155, 160, "Trans", '#833792')
    ]
    for start, end, label, color in extra_bars:
        fig.add_shape(type="rect", x0=start - 0.5, x1=end - 0.5,
                      y0=y_max * 1.08, y1=y_max * 1.16,   # <-- thicker
                      line=dict(width=0), fillcolor=color, opacity=1)
        fig.add_annotation(x=(start + end) / 2 - 0.5, y=(y_max * 1.08 + y_max * 1.15) / 2,
                           text=label, showarrow=False, font=dict(size=12, color="white"))

    def add_horizontal_bar(start, end, y, label, color, opacity=1):
        fig.add_shape(type="rect", x0=start - 0.5, x1=end - 0.5,
                      y0=y, y1=y * (1.22 / 1.15),   # slightly thicker too
                      line=dict(width=0), fillcolor=color, opacity=opacity)
        fig.add_annotation(x=(start + end) / 2 - 0.5, y=(y + y * (1.18 / 1.11)) / 2,
                           text=label, showarrow=False, font=dict(size=12, color="white"))

    horizontal_bars = [
        (0, 80, y_max * 1.15, "Clustered", 'black'),
        (80, 160, y_max * 1.15, "Non-Clustered", 'lightgrey')
    ]
    for bar in horizontal_bars:
        add_horizontal_bar(*bar)

    tickvals = np.arange(161)
    section_1 = ['1-10Kb very early', '1-10Kb early', '1-10Kb medium', '1-10Kb late', '1-10Kb very late']
    section_2 = ['10-100Kb very early', '10-100Kb early', '10-100Kb medium', '10-100Kb late', '10-100Kb very late']
    section_3 = ['100Kb-1Mb very early', '100Kb-1Mb early', '100Kb-1Mb medium', '100Kb-1Mb late', '100Kb-1Mb very late']
    section_4 = ['1Mb-10Mb very early', '1Mb-10Mb early', '1Mb-10Mb medium', '1Mb-10Mb late', '1Mb-10Mb very late']
    section_5 = ['>10Mb very early', '>10Mb early', '>10Mb medium', '>10Mb late', '>10Mb very late']
    repeated_list = (section_1 + section_2 + section_3 + section_4 + section_5) * 3 + \
                    ['very early', 'early', 'medium', 'late', 'very late'] + \
                    (section_1 + section_2 + section_3 + section_4 + section_5) * 3 + \
                    ['very early', 'early', 'medium', 'late', 'very late']

    fig.update_layout(
        xaxis=dict(tickmode='array', tickvals=tickvals, ticktext=repeated_list, tickangle=-90),
        height=400, width=1100,
        yaxis_title="Fraction of Total Mutations",
        xaxis_title="Mutation Types",
        xaxis_tickangle=-90,
        yaxis=dict(range=[0, y_max * 1.22]),   # raise range slightly
        showlegend=False,
        paper_bgcolor='white',
        plot_bgcolor='white'
    )

    fig.update_layout(
        font=dict(
            family="Arial",
            size=9,
            color="black"
        )
    )

    if show:
        fig.show()
        fig.write_image(f"{output_dir}/{title}.jpg", format='jpg', scale=3)





def plot_rt_sig_2d(mut_spectrum, output_dir, sig):
    mut_spectrum_reindexed = mut_spectrum.copy()
    mut_spectrum_reindexed.index = mut_spectrum_reindexed['MutationType']
    mut_spectrum_reindexed['rt'] = mut_spectrum_reindexed['MutationType'].str.split('_').str[-1]
    mut_spectrum_reindexed['class'] = (
        mut_spectrum_reindexed['MutationType'].str.split('_').str[1] + '_' +
        mut_spectrum_reindexed['MutationType'].str.split('_').str[2]
    )
    mut_spectrum_reindexed['class'] = mut_spectrum_reindexed['class'].str.replace(
    r'_(very early|early|medium|late|very late)$', '', regex=True)
    mut_spectrum_reindexed['is_clustered'] = mut_spectrum_reindexed['MutationType'].str.split('_').str[0]=='clustered'
    mut_spectrum_reindexed['col'] = [colDict[prefix] for prefix in mut_spectrum_reindexed['MutationType'].str.split('_').str[1]]
    mut_spectrum_reindexed_nonclustered = mut_spectrum_reindexed[mut_spectrum_reindexed['is_clustered']==False]
    plt.figure(figsize=(15, 4))
    # Convert classes to integer positions
    x = pd.factorize(mut_spectrum_reindexed_nonclustered['class'])[0]
    y = mut_spectrum_reindexed_nonclustered['rt'].values
    s = mut_spectrum_reindexed_nonclustered['Total Count'].values * 1000
    ax = plt.gca()
    # Make sure layout is known before coordinate transform
    plt.scatter(x, y, s=s, alpha=0)  # invisible points to set axis scale
    plt.gcf().canvas.draw()
    # Convert marker size from points to data units to compute half-width
    side_pts = np.sqrt(s)  # full side length in points
    trans = ax.transData.inverted()
    half_widths_data = []
    for xi, w_pts in zip(x, side_pts):
        x0_disp = ax.transData.transform((xi, 0))
        x1_disp = (x0_disp[0] + w_pts / 2, x0_disp[1])  # 
        x1_data = trans.transform(x1_disp)
        half_widths_data.append(x1_data[0] - xi)
    half_widths_data = np.array(half_widths_data)
    # Shift x to the RIGHT by half width → left edge aligns at original x
    x_shifted = x + half_widths_data
    # Actual scatter plot
    plt.scatter(
    x=x_shifted,
    y=y,
    s=s,
    alpha=1,
    c=mut_spectrum_reindexed['col'],
    edgecolors='none',
    marker='s'
    )
    # Axis formatting
    # Axis formatting
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    
    # Remove x-axis labels and ticks completely
    ax.set_xticks(x)                      # one tick per class
    ax.set_xticklabels([])                # hide labels, keep ticks for grid
    
    #Enlarge y-axis label and title 2.5×
    base_fontsize = 12
    ax.set_ylabel('RFD', fontsize=base_fontsize * 2.5)
    ax.set_title(sig, fontsize=base_fontsize * 2.5)
    
    # Also enlarge y-tick labels
    ax.tick_params(axis='y', labelsize=base_fontsize * 2.5)
    
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
    output_path = f"{output_dir}/{sig}_scatter_all_leftaligned.pdf"
    plt.savefig(output_path, format='pdf', dpi=300)
