# Load data -------------------------------------------------------------------#
# ps_rel <- readRDS("RDS/02_ps_rel_bac_norm.rds")

# PCoA with weighted UniFrac Analysis -----------------------------------------#
pcoa_plot <- ps_pcoa(ps = ps_rel, dist.metric = "wunifrac", col_name = data_00$col_name)

ggsave("04_pcoa_test.png",
       plot = pcoa_plot,
       width = 15,
       height = 5)