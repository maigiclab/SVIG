import pandas as pd
import numpy as np
import plotly.express as px
import matplotlib.pyplot as plt
from matplotlib import gridspec

order = ["L_L", "L_R", "R_L", "R_R"]

color_map = {
    # Deletions - shades of red
    'Clustered Deletion (L_L)': '#c61033', # Light red
    'Clustered Deletion (L_R)': '#c61033',  # Medium red
    'Clustered Deletion (R_L)': '#c61033',  # Dark red
    'Clustered Deletion (R_R)': '#c61033', # Bright red
    'Non-clustered Deletion (L_L)': '#c61033',  # Light red
    'Non-clustered Deletion (L_R)': '#c61033',  # Medium red
    'Non-clustered Deletion (R_L)': '#c61033',  # Dark red
    'Non-clustered Deletion (R_R)': '#c61033',  # Bright red
    
    # Inversions - shades of blue
    'Clustered Inversion (L_L)': '#3b74ac',  # Light blue
    'Clustered Inversion (L_R)': '#3b74ac',  # Medium blue
    'Clustered Inversion (R_L)': '#3b74ac',  # Dark blue
    'Clustered Inversion (R_R)': '#3b74ac',  # Bright blue
    'Non-clustered Inversion (L_L)': '#3b74ac',  # Light blue
    'Non-clustered Inversion (L_R)': '#3b74ac',  # Medium blue
    'Non-clustered Inversion (R_L)': '#3b74ac',  # Dark blue
    'Non-clustered Inversion (R_R)': '#3b74ac',  # Bright blue
    
    # Tandem Duplications - shades of green
    'Clustered Tandem Duplication (L_L)': '#44a23c',  # Light green
    'Clustered Tandem Duplication (L_R)': '#44a23c',  # Medium green
    'Clustered Tandem Duplication (R_L)': '#44a23c',  # Dark green
    'Clustered Tandem Duplication (R_R)': '#44a23c',  # Bright green
    'Non-clustered Tandem Duplication (L_L)': '#44a23c',  # Light green
    'Non-clustered Tandem Duplication (L_R)': '#44a23c',  # Medium green
    'Non-clustered Tandem Duplication (R_L)': '#44a23c',  # Dark green
    'Non-clustered Tandem Duplication (R_R)': '#44a23c',  # Bright green
    
    # Translocations - shades of purple
    'Clustered Translocation (L_L)': '#833792', # Light purple
    'Clustered Translocation (L_R)': '#833792',  # Medium purple
    'Clustered Translocation (R_L)': '#833792',  # Dark purple
    'Clustered Translocation (R_R)': '#833792',  # Bright purple
    'Non-clustered Translocation (L_L)': '#833792',  # Light purple
    'Non-clustered Translocation (L_R)': '#833792',  # Medium purple
    'Non-clustered Translocation (R_L)': '#833792',  # Dark purple
    'Non-clustered Translocation (R_R)': '#833792',  # Bright purple
}

colDict = {}
colDict['del'] = '#c61033'
colDict['tds'] = '#44a23c'
colDict['inv'] = '#3b74ac'
colDict['trans'] = '#833792'

# barplot plotting of the full RFD signature
# ----------------- Plotting function -----------------
def plot_mut_spectrum(mutation_vector, sorted_array, color_map, title, output_dir, show=True, showPdf=False):
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
        subcat = "_".join(mutation.split("_")[-2:])
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
    for i in range(0, 129, 4):
        fig.add_vline(x=i - 0.5, line_dash="dot", line_color="grey", opacity=0.6, line_width=0.3)

    extra_bars = [
        (0, 20, "Deletions", '#c61033'),
        (20, 40, "Duplications", '#44a23c'),
        (40, 60, "Inversions", '#3b74ac'),
        (60, 64, "Trans", '#833792'),
        (64, 84, "Deletions", '#c61033'),
        (84, 104, "Duplications", '#44a23c'),
        (104, 124, "Inversions", '#3b74ac'),
        (124, 128, "Trans", '#833792')
    ]

    for start, end, label, color in extra_bars:
        fig.add_shape(type="rect", x0=start - 0.5, x1=end - 0.5,
                      y0=y_max * 1.10, y1=y_max * 1.15,
                      line=dict(width=0), fillcolor=color, opacity=1)
        fig.add_annotation(x=(start + end) / 2 - 0.5, y=(y_max * 1.10 + y_max * 1.15) / 2,
                           text=label, showarrow=False, font=dict(size=10, color="white"))

    def add_horizontal_bar(start, end, y, label, color, opacity=1):
        fig.add_shape(type="rect", x0=start - 0.5, x1=end - 0.5,
                      y0=y, y1=y * (1.20 / 1.15),
                      line=dict(width=0), fillcolor=color, opacity=opacity)
        fig.add_annotation(x=(start + end) / 2 - 0.5, y=(y + y * (1.15 / 1.10)) / 2,
                           text=label, showarrow=False, font=dict(size=10, color="white"))

    horizontal_bars = [
        (0, 64, y_max * 1.15, "Clustered", 'black'),
        (64, 128, y_max * 1.15, "Non-Clustered", 'lightgrey')
    ]
    for bar in horizontal_bars:
        add_horizontal_bar(*bar)

    fig.update_layout(
        yaxis_title="Fraction of Total Mutations",
        xaxis_title="Mutation Types",
        xaxis_tickangle=-90,
        yaxis=dict(range=[0, y_max * 1.20]),
        showlegend=False
    )

    section_labels = ['1-10Kb', '10-100Kb', '100Kb-1Mb', '1Mb-10Mb', '>10Mb']
    orientations = ['L_L', 'L_R', 'R_L', 'R_R']
    combined_sections = [f"{sec} {ori}" for sec in section_labels for ori in orientations]
    repeated_labels = combined_sections * 3 + orientations + combined_sections * 3 + orientations

    fig.update_layout(
        xaxis=dict(tickmode='array', tickvals=np.arange(128), ticktext=repeated_labels, tickangle=-90),
        height=400, width=1100
    )

    fig.update_layout(
        showlegend=False,
        paper_bgcolor='white',  # Set background color to white
        plot_bgcolor='white'    # Set plot area background color to white
    )

    fig.update_layout(
        font=dict(
            family="Arial",
            size=9,
            color="black"
        )
    )

    if show:
        if showPdf:
            fig.show()
            print('writing PDF')
            pdf_path = f"{output_dir}/{title}.pdf"
            print(pdf_path)
            fig.write_image(pdf_path, format='pdf', width=1100, height=400, scale=3)
        else: 
            fig.show()
            print(f"{output_dir}/{title}.jpg")
            fig.write_image(f"{output_dir}/{title}.jpg", format='jpg', scale=3)


# like above, but has the additional functionality to show a subset of channels, eg only_nonclustered_td
def plot_mut_spectrum2(
    mutation_vector, sorted_array, color_map, title, output_dir,
    show=True, showPdf=False, only_nonclustered_td=False
):
    import pandas as pd
    import numpy as np
    import plotly.express as px
    import plotly.graph_objects as go

    # -----------------------------------------------------------
    # Helper: classify mutation
    # -----------------------------------------------------------
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
        subcat = "_".join(mutation.split("_")[-2:])
        return f"{cluster} {mtype} ({subcat})"

    df = mutation_vector.copy()

    # -----------------------------------------------------------
    # FILTER MODE: Only non-clustered tandem duplications
    # -----------------------------------------------------------
    if only_nonclustered_td:

        df = df[
            df["MutationType"].str.startswith("non-clustered_") &
            df["MutationType"].str.contains("_tds_") &
            ~df["MutationType"].str.contains('tds_>10Mb')
        ]
        df["Group"] = df["MutationType"].apply(categorize_mutation)

        # Update category list to match filtered data
        sorted_array = [c for c in sorted_array if c in df["MutationType"].unique()]

        # Extract prefix (SV size bin) and suffix (orientation)
        df["prefix"] = df["MutationType"].str.extract(r"tds_(.+?)_[LR]_[LR]")
        df["suffix"] = df["MutationType"].str.extract(r"(_[LR]_[LR])$")

        # Sort by prefix first, then suffix in L_L, L_R, R_L, R_R order
        suffix_order = ["L-L", "L-R", "R-L", "R-R"]
        df["suffix"] = pd.Categorical(
            df["suffix"].str[1:].str.replace("_", "-", regex=False),
            categories=suffix_order,
            ordered=True
        )
        
        df.sort_values(["prefix", "suffix"], inplace=True)

        # Reset index so that x positions are clean 0..N-1
        df.reset_index(drop=True, inplace=True)

        # -----------------------------------------------------------
        # TD-only plot
        # -----------------------------------------------------------
        fig = px.bar(
            df,
            x=df.index,  # x = numeric positions
            y="Total Count",
            color="Group",
            labels={"Total Count": "Mutation Count"},
            color_discrete_map=color_map
        )

        fig.update_traces(marker_line_width=0)
        fig.update_layout(title_x=0.5)

        y_max = df["Total Count"].max() if len(df) > 0 else 1

        # -----------------------------------------------------------
        # Vertical dashed lines between categories
        # -----------------------------------------------------------
        for i in range(0, len(df), 4):
            fig.add_vline(
                x=i - 0.5,
                line_dash="dot",
                line_color="grey",
                opacity=0.6,
                line_width=0.6
            )

        # -----------------------------------------------------------
        # Replace x-ticks with short suffix labels
        # -----------------------------------------------------------
        fig.update_xaxes(
            tickmode="array",
            tickvals=list(df.index),
            ticktext=df["suffix"].astype(str).tolist(),
            tickangle=-45
        )

        # -----------------------------------------------------------
        # Add prefix labels above groups of 4 bars
        # -----------------------------------------------------------
        for prefix in df["prefix"].unique():
            group_idx = df.index[df["prefix"] == prefix].tolist()
            if len(group_idx) == 0:
                continue
            start = min(group_idx)
            end = max(group_idx)
            mid = (start + end) / 2

            fig.add_annotation(
                x=mid,
                y=y_max * 1.05,
                text=prefix,
                showarrow=False,
                font=dict(size=11, color="black")
            )

        fig.update_layout(
            height=450,
            width=739,
            yaxis=dict(range=[0, y_max * 1.20]),
            xaxis_title="",
            yaxis_title="Mutation Count",
            plot_bgcolor="white",
            paper_bgcolor="white",
            showlegend=False,
            font=dict(family="Arial", size=10, color="black")
        )
        fig.update_layout(
            margin=dict(l=0, r=0, t=10, b=0),
        )

    # =======================================================================
    # FULL COMPLEX MODE (original 128-category mutation spectrum)
    # =======================================================================
    else:
        df["Group"] = df["MutationType"].apply(categorize_mutation)
        df["Color"] = df["Group"].map(color_map)
        df["MutationType"] = pd.Categorical(df["MutationType"],
                                            categories=sorted_array,
                                            ordered=True)
        df.sort_values("MutationType", inplace=True)
        df.reset_index(drop=True, inplace=True)

        fig = px.bar(
            df,
            x="MutationType",
            y="Total Count",
            color="Group",
            title=title,
            labels={"Total Count": "Mutation Count", "MutationType": "Mutation Type"},
            color_discrete_map=color_map,
            text_auto=False,
            category_orders={"MutationType": sorted_array}
        )

        fig.update_traces(marker_line_width=0)
        fig.update_layout(title_x=0.5)

        y_max = df["Total Count"].max()

        # Vertical dashed lines
        for i in range(0, 129, 4):
            fig.add_vline(x=i - 0.5, line_dash="dot",
                          line_color="grey", opacity=0.6, line_width=0.3)

        # Header rectangles + labels (unchanged from your original)
        extra_bars = [
            (0, 20, "Deletions", '#c61033'),
            (20, 40, "Duplications", '#44a23c'),
            (40, 60, "Inversions", '#3b74ac'),
            (60, 64, "Trans", '#833792'),
            (64, 84, "Deletions", '#c61033'),
            (84, 104, "Duplications", '#44a23c'),
            (104, 124, "Inversions", '#3b74ac'),
            (124, 128, "Trans", '#833792')
        ]

        for start, end, label, color in extra_bars:
            fig.add_shape(
                type="rect",
                x0=start - 0.5, x1=end - 0.5,
                y0=y_max * 1.10, y1=y_max * 1.15,
                line=dict(width=0),
                fillcolor=color, opacity=1
            )
            fig.add_annotation(
                x=(start + end) / 2 - 0.5,
                y=(y_max * 1.10 + y_max * 1.15) / 2,
                text=label,
                showarrow=False,
                font=dict(size=10, color="white")
            )

        # Clustered / Non-clustered grey/black bars
        def add_horizontal_bar(start, end, y, label, color, opacity=1):
            fig.add_shape(type="rect",
                          x0=start - 0.5, x1=end - 0.5,
                          y0=y, y1=y * (1.20 / 1.15),
                          line=dict(width=0), fillcolor=color, opacity=opacity)
            fig.add_annotation(
                x=(start + end) / 2 - 0.5,
                y=(y + y * (1.15 / 1.10)) / 2,
                text=label, showarrow=False,
                font=dict(size=10, color="white")
            )

        horizontal_bars = [
            (0, 64, y_max * 1.15, "Clustered", 'black'),
            (64, 128, y_max * 1.15, "Non-Clustered", 'lightgrey')
        ]
        for bar in horizontal_bars:
            add_horizontal_bar(*bar)

        # Custom tick labels for full spectrum
        section_labels = ['1-10Kb', '10-100Kb', '100Kb-1Mb', '1Mb-10Mb', '>10Mb']
        orientations = ['L_L', 'L_R', 'R_L', 'R_R']
        combined_sections = [f"{sec} {ori}" for sec in section_labels for ori in orientations]
        repeated_labels = combined_sections * 3 + orientations + combined_sections * 3 + orientations

        fig.update_layout(
            xaxis=dict(
                tickmode='array',
                tickvals=np.arange(128),
                ticktext=repeated_labels,
                tickangle=-90
            ),
            height=400,
            width=1100,
            yaxis=dict(range=[0, y_max * 1.20]),
            showlegend=False,
            plot_bgcolor='white',
            paper_bgcolor='white',
            font=dict(family="Arial", size=9, color="black")
        )

    # -----------------------------------------------------------
    # SAVE / SHOW
    # -----------------------------------------------------------
    if show:
        if showPdf:
            fig.show()
            pdf_path = f"{output_dir}/{title}.pdf"
            print("writing PDF:", pdf_path)
            if only_nonclustered_td:
                fig.write_image(pdf_path, format='pdf', width=655, height=400, scale=3)
            else:
                fig.write_image(pdf_path, format='pdf', width=1100, height=400, scale=3)
        else:
            fig.show()
            jpg_path = f"{output_dir}/{title}.jpg"
            print(jpg_path)
            fig.write_image(jpg_path, format='jpg', scale=3)

# function for drawing RFD signatures in 2 dimensions
def plot_square_small(mut_spectrum_reindexed_sel, output_dir_corr, sig):
    # Example square plot
        plt.figure(figsize=(5, 4))
    
        mut_spectrum_reindexed_sel["rfd"] = pd.Categorical(
            mut_spectrum_reindexed_sel["rfd"],
            categories=order,
            ordered=True
        )
    
        # Convert classes to integer positions
        x = pd.factorize(mut_spectrum_reindexed_sel['class'])[0]
        y = mut_spectrum_reindexed_sel['rfd'].cat.codes
        s = mut_spectrum_reindexed_sel['Total Count'].values * 1000
        
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
            x1_disp = (x0_disp[0] + w_pts / 2, x0_disp[1])  # 👉 half width!
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
            c=mut_spectrum_reindexed_sel['col'],
            edgecolors='none',
            marker='s'
        )
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels(order)
    
        # Axis formatting
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
        ax.set_xticks(x)
        ax.set_xticklabels(mut_spectrum_reindexed_sel['class'], rotation=45, ha='right')
        
        plt.ylabel('', fontsize=12)
        plt.title(sig, fontsize=14)
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        output_path = f"{output_dir_corr}/{sig}_scatter_leftaligned.pdf"
        plt.savefig(output_path, format='pdf', dpi=300)
        plt.show()

def plot_square_small2(
    mut_spectrum_reindexed_sel,
    output_dir_corr,
    sig,
    add_row_boxplot=False,
    barplot_col='gray',
    disable_y_tick_labels=False
):
    # --- Figure & layout ---
    if add_row_boxplot:
        fig = plt.figure(figsize=(6.5, 4))
        gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1], wspace=0.05)
        ax = fig.add_subplot(gs[0])
        ax_box = fig.add_subplot(gs[1], sharey=ax)
    else:
        fig, ax = plt.subplots(figsize=(5, 4))

    # --- Category ordering ---
    mut_spectrum_reindexed_sel["rfd"] = pd.Categorical(
        mut_spectrum_reindexed_sel["rfd"],
        categories=order,
        ordered=True
    )

    # --- Coordinates & sizes ---
    x = pd.factorize(mut_spectrum_reindexed_sel["class"])[0]
    y = mut_spectrum_reindexed_sel["rfd"].cat.codes
    s = mut_spectrum_reindexed_sel["Total Count"].values * 1000

    # Invisible scatter to define axis scale
    ax.scatter(x, y, s=s, alpha=0)
    fig.canvas.draw()

    # Convert marker size (points) → data units
    side_pts = np.sqrt(s)
    trans = ax.transData.inverted()
    half_widths_data = []

    for xi, w_pts in zip(x, side_pts):
        x0_disp = ax.transData.transform((xi, 0))
        x1_disp = (x0_disp[0] + w_pts / 2, x0_disp[1])
        x1_data = trans.transform(x1_disp)
        half_widths_data.append(x1_data[0] - xi)

    half_widths_data = np.array(half_widths_data)
    x_shifted = x + half_widths_data

    # --- Main square plot ---
    ax.scatter(
        x=x_shifted,
        y=y,
        s=s,
        c=mut_spectrum_reindexed_sel["col"],
        edgecolors="none",
        marker="s"
    )

    ax.set_yticks(range(len(order)))
    if disable_y_tick_labels:
        ax.set_yticklabels([])
    else:
        ax.set_yticklabels(order)
        
    x_unique = np.unique(x)
    class_unique = mut_spectrum_reindexed_sel["class"].iloc[
        np.unique(x, return_index=True)[1]
    ]

    ax.set_xticks(x)                      # one tick per class
    ax.set_xticklabels([])                # hide labels, keep ticks for grid
    
    ax.grid(True, which='both', axis='both', alpha=0.3)




    ax.set_ylabel("", fontsize=12)
    ax.set_title(sig, fontsize=14)
    ax.grid(True, alpha=0.3)

    # --- Optional row-sum boxplot ---
    if add_row_boxplot:
        row_sums = (
            mut_spectrum_reindexed_sel
            .groupby("rfd", observed=True)["Total Count"]
            .sum()
            .reindex(order)
        )

        ax_box.barh(
            y=np.arange(len(order)),
            width=row_sums.values,
            height=0.3,
            color=barplot_col
        )

        ax_box.set_xlabel("Weight (norm)", fontsize=10)
        ax_box.grid(True, axis="x", alpha=0.3)

        # Clean alignment with main axis
        ax_box.tick_params(axis="y", left=False, labelleft=False)
        ax_box.invert_yaxis()  # match categorical order
        ax_box.spines["top"].set_visible(False)
        ax_box.spines["right"].set_visible(False)
        ax_box.spines["left"].set_visible(False)

    # --- Output ---
    plt.tight_layout()
    output_path = f"{output_dir_corr}/{sig}_scatter_leftaligned.pdf"
    plt.savefig(output_path, format="pdf", dpi=300)
    plt.show()

def plot_square_big(mut_spectrum_reindexed, output_dir_corr, sig):
    plt.figure(figsize=(12, 4))
    
    # Convert classes to integer positions
    x = pd.factorize(mut_spectrum_reindexed['class'])[0]
    y = mut_spectrum_reindexed['rfd'].cat.codes
    s = mut_spectrum_reindexed['Total Count'].values * 1000
    
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
        x1_disp = (x0_disp[0] + w_pts / 2, x0_disp[1])  # 👉 half width!
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
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order)
    
    # Axis formatting
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    
    #  Remove x-axis labels and ticks completely
    ax.set_xticks(x)                      # one tick per class
    ax.set_xticklabels([])                # hide labels, keep ticks for grid
    
    ax.grid(True, which='both', axis='both', alpha=0.3)
    
    # Enlarge y-axis label and title 2.5×
    base_fontsize = 12
    ax.set_ylabel('RFD', fontsize=base_fontsize * 2.5)
    ax.set_title(sig, fontsize=base_fontsize * 2.5)
    
    # Also enlarge y-tick labels
    ax.tick_params(axis='y', labelsize=base_fontsize * 2.5)
    
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    plt.tight_layout()
    output_path = f"{output_dir_corr}/{sig}_scatter_all_leftaligned.pdf"
    plt.savefig(output_path, format='pdf', dpi=300)
    plt.show()


def heatmap_style(g):
    """
    Apply consistent styling to a seaborn clustermap,
    reposition colorbar, make it smaller, and add vertical lines every 16 columns.
    """
    if hasattr(g, 'ax_heatmap'):  # clustermap
        ax = g.ax_heatmap
        cbar = g.cax

        # Shrink heatmap slightly to leave space for colorbar
        pos = ax.get_position()
        ax.set_position([pos.x0, pos.y0, pos.width*0.85, pos.height])

        # Reposition colorbar (smaller width)
        cbar.set_position([0.88, pos.y0, 0.015, pos.height])
        cbar.tick_params(width=0.5, length=4, labelsize=10)

        # Ticks styling
        ax.tick_params(axis='both', width=0.5, length=4)

        # --- Add vertical lines every 16 columns ---
        n_cols = len(ax.get_xticklabels())
        #for x in range(4, n_cols, 4):
        #    ax.axvline(x, color='black', linewidth=1)
    else:
        ax = g
        cbar = ax.collections[0].colorbar
        cbar.ax.set_position([0.88, 0.2, 0.015, 0.6])
        ax.tick_params(axis='both', width=0.5, length=4)
        cbar.ax.tick_params(width=0.5, length=4, labelsize=10)

        # vertical lines
        n_cols = len(ax.get_xticklabels())
        for x in range(4, n_cols, 4):
            ax.axvline(x, color='black', linewidth=1)
