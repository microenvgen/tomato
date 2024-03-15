

# =====LIBRARIES================================================================
r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)

library(phyloseq)

load("cdpcoa_input.RData")
                    
ps.dpcoa <- ordinate(ps_si3, "DPCoA")          
              
save.image("cdpcoa_results.RData")





