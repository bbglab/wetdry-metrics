#!/usr/local/bin/python

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import statsmodels.api as sm
import numpy as np
from datetime import datetime

from matplotlib.backends.backend_pdf import PdfPages
from utils_definitions import *


cols_dry_only = [
    'DryLab>>On-target unique molecules percentage',
    'DryLab>>SSC to DSC',
    'DryLab>>Raw to DSC',
    'DryLab>>Amplification bias on/off',
]

cols_complete = cols_dry_only + [
    'Combined>>Unique molecules sequenced vs qPCR',
    'Combined>>Recovery input to depth',
    'WetLab>>Recovery Input to ligation',
    'WetLab>>Recovery Input to ligation qPCR'
]


def plot_metrics(metrics_df, output_plot_dir) : 

    ### Plot predictive metrics per cohort ###

    filename = f'{output_plot_dir}/metrics_plots.pdf'

    with PdfPages(filename) as pdf:

        if 'Combined>>Unique molecules sequenced vs qPCR' in metrics_df.columns :
            columns_2plot = cols_complete
            by = 'BATCH'

        else :
            columns_2plot = cols_dry_only
            by = 'Run ID'

        for var in columns_2plot:
            complete_to_plot = metrics_df

            # Plot per project
            plt.figure(figsize=(10, 6))
            sns.boxplot(data=complete_to_plot,
                        x=by,
                        y=var,
                        showfliers = False,
                        legend = False,
                    )
            sns.stripplot(data=complete_to_plot,
                        x=by,
                        y=var,
                        jitter=True,
                        color='black',
                        size=3,
                        legend = False,
                        )
            plt.title(var)
            plt.xticks(rotation = 90)
            # plt.legend(bbox_to_anchor = (1.1, 1))
            plt.tight_layout()
            
            # Save the plot to the PDF
            pdf.savefig()
            #plt.show()
            plt.close()
                    
        ### Final message ###

    print('Plotting done')