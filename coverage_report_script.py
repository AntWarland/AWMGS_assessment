import pandas as pd

# set global variables #
patient_id = "A"
depth_col = "Total_Depth"
locus_col = "Locus"
region_col = "region"

# read in data #
names=["chr", "start", "end", region_col]
df_roi = pd.read_csv("lung.bed", names=names, delimiter="\t")
df_depth = pd.read_csv("depth_of_coverage.csv")

# split piped columns into separate columns #
df_roi[["gene", "code", "exon"]] = df_roi[region_col].str.split("|", n=2, expand=True)
df_roi[region_col] = df_roi[region_col].str.replace("|", "_")
df_depth[["chr", "loc"]] = df_depth[locus_col].str.split(":", n=1, expand=True)

# convert location strings to integers #
df_depth["location"] = pd.to_numeric(df_depth["loc"], errors="coerce", downcast="integer")

# drop any depth NAs #
df_depth.dropna(subset=[depth_col], inplace=True)

def calc_depth_stats(df, depth_column, name):
    """
    Calculates mean coverage and percentage of bases with coverage > 250 for a given df

    Args:
        df (df): dataframe with coverage depth values
        depth_column (str): depth column name
        region name (str): slice/region/exon description
    
    Returns:
        dataframe containing region name, average coverage and percentage
    """
    # initialise summary stats df #
    stats_df = pd.DataFrame(columns=["region", "average_coverage", "%250x"])
    
    # calculate mean coverage #
    depth_series = df[depth_column]
    mean = int(round(depth_series.mean()))
    
    # calculate % of bases with coverage > 250x #
    length = len(depth_series)
    length_250 = len(depth_series[depth_series>=250])
    percent_250 = round(100*(length_250/length),1)
    
    # store values in stats df # 
    stats_df.loc[len(df)] = [name, mean, percent_250]
    return stats_df

def summarise_depths(rois, depth):
    """
    Creates a summary table of per-base depth statistics for regions of interest
    
    Args: 
        roi: dataframe with regions of interest, including locus start and end positions
        depth: dataframe with per base depth of coverage for regions of interest

    Returns:
        Summary table containing coverage and % of bases confidently called over all regions,
        gene averages and specific regions
    """
    # initialise empty dfs # 
    summary_df = pd.DataFrame()
    full_df = pd.DataFrame()
    
    # aggregate per region stats in summary_df #
    for i, row in rois.iterrows():
        chr_slice = depth.loc[depth["chr"]==row["chr"], :].copy()
        locus_slice = chr_slice.loc[(chr_slice["location"]>=row["start"]) 
                                    & (chr_slice["location"]<=row["end"]), :].copy()
        locus_stats = calc_depth_stats(locus_slice, depth_col, row[region_col])
        locus_slice["gene_code"] = f"{row.gene}_{row.code}"
        full_df = pd.concat([full_df, locus_slice], ignore_index=True)
        summary_df = pd.concat([summary_df, locus_stats], ignore_index=True)
    
    # list genes_codes #
    gene_codes = list(full_df["gene_code"].unique())
    
    # add per-gene stats to summary table #
    for i in gene_codes:
        gene_code_slice = full_df.loc[full_df["gene_code"]==i].copy()
        name = i + "_all"
        gene_code_stats = calc_depth_stats(gene_code_slice, depth_col, name)
        summary_df = pd.concat([summary_df, gene_code_stats], ignore_index=True)
    
    # sort summary table #
    summary_df.sort_values(["region"], inplace=True)
    
    # get whole panel stats and add to top of summary table #
    all_stats = calc_depth_stats(full_df, depth_col, "whole_panel")
    summary_df = pd.concat([all_stats, summary_df], ignore_index=True)

    return summary_df

# execute summary report function #
coverage = summarise_depths(df_roi, df_depth)

# export dataframe to tsv file # 
coverage.to_csv(f"patient_{patient_id}_coverage_report.tsv", sep='\t', index=False, header=True)