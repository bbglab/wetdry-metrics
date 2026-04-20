#!/usr/local/bin/python

target_amplification = 18 # TWS
size_whole_genome = 3088286376
adapter_size = 150
read_size = 150
umi_size = 9
adapter_index_size = 180
genome_weight = 0.003 # ng
genome_size = 3e9
avogadro = 6.02e23
avogadro_in_f_scale = avogadro / 1e15
bp_weight_in_ng = genome_weight / genome_size
max_contribution_per_read = 300 - 18

informative_columns_drylab = [
                            "SAMPLE_ID",
                            "Run ID", 
                            "PATH",
                            
                            "DryLab>>GBs analysed",
                            "DryLab>>Total GBs for optimal",
                            'DryLab>>GBs (analyzed-optimal)',


                            "DryLab>>Depth", 
                            
                            'DryLab>>On-target unique molecules percentage',
                            
                            "DryLab>>Percentage of duplicates on-target", 
                            "DryLab>>Family size",
                            "DryLab>>SSC to DSC",
                            "DryLab>>Raw to DSC",
                            "DryLab>>Amplification bias on/off",
                            'FamMetrics>>on_target.duplex.raw_x_dscs',
                            'FamMetrics>>off_target.duplex.raw_x_dscs',
                            
                            'FamMetrics>>on_target.raw_reads',
                            'FamMetrics>>off_target.raw_reads',
                            ]

informative_columns = [
                        "DryLab ID", 
                        "WetLab ID",
                        "BATCH", 
                        "Run ID", 
                        "PATH",
                        "WetLab>>Panel",
                        "WetLab>>Number of captures",
                        
                        "WetLab>>Input (ng)",
                        "WetLab>>fmol to PCR1",
                        'WetLab>>Recovery Input to ligation',
                        'WetLab>>Recovery Input to ligation qPCR',
                        
                        'WetLab>>qPCR unique molecules',

                        "DryLab>>GBs analysed",
                        "DryLab>>Total GBs for optimal",
                        'DryLab>>GBs (analyzed-optimal)',

                        "DryLab>>Depth", # depth
                        
                        "Combined>>Unique molecules sequenced vs qPCR",
                        "Combined>>Recovery input to depth",

                        'DryLab>>On-target unique molecules percentage',

                        "DryLab>>Percentage of duplicates on-target",
                        "DryLab>>Family size", # family size
                        "DryLab>>SSC to DSC",
                        "DryLab>>Raw to DSC",
                        "DryLab>>Amplification bias on/off",

                        'FamMetrics>>on_target.duplex.duplex_raw_x_dscs',
                        'FamMetrics>>off_target.duplex.raw_x_dscs',
                        
                        'FamMetrics>>on_target.raw_reads',
                        'FamMetrics>>off_target.raw_reads',
                        
                        ]   
