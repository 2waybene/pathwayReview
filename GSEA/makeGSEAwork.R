sample.ds = system.file('extdata', 'Leukemia_hgu95av2.gct', package = 'GSEA', mustWork = TRUE)


gseares <- GSEA(input.ds = "GSEA/example/Leukemia_hgu95av2.gct",
                input.cls = "GSEA/example/Leukemia.cls",
                input.chip = "GSEA/example/Human_AFFY_HG_U95_MSigDB_7_0_final.chip",
                gs.db = "GSEA/GSEA_db/c2.cp.kegg.v7.0.symbols.gmt",
                
                #gs.db = "GSEA/GSEA_db/c2.cp.kegg.v7.0.entrez.gmt",
                #output.directory = "gsea_stored_results/",
                
                output.directory = "GSEA/",
                doc.string= "kegg",
                #   non.interactive.run   = F,               
                reshuffling.type      = "sample.labels", 
                nperm                 = 1000,
                weighted.score.type   =  1,            
                nom.p.val.threshold   = -1,            
                fwer.p.val.threshold  = -1,            
                fdr.q.val.threshold   = 0.25,          
                topgs                 = 20,            
                adjust.FDR.q.val      = F,             
                gs.size.threshold.min = 1,             
                gs.size.threshold.max = 500,           
                reverse.sign          = F,             
                preproc.type          = 0,             
                random.seed           = 3338,          
                perm.type             = 0,             
                fraction              = 1.0,           
                replace               = F,             
                save.intermediate.results = F,         
                #   OLD.GSEA              = F,             
                use.fast.enrichment.routine = T,
                
                collapse.dataset = TRUE, collapse.mode = 'max')
                
              #  )


gseaEG <-GSEA(
  
#     input.ds = system.file('extdata', 'Leukemia_hgu95av2.gct', package = 'GSEA', mustWork = TRUE),
#     input.cls = system.file('extdata', 'Leukemia.cls', package = 'GSEA', mustWork = TRUE),
#     input.chip = system.file('extdata', 'Human_AFFY_HG_U95_MSigDB_7_0_final.chip',
#                              package = 'GSEA', mustWork = TRUE), 
     
     input.ds = "GSEA/example/Leukemia_hgu95av2.gct",
     input.cls = "GSEA/example/Leukemia.cls",
     input.chip = "GSEA/example/Human_AFFY_HG_U95_MSigDB_7_0_final.chip",
     
     
     #gs.db = system.file('extdata','h.all.v7.0.symbols.gmt', package = 'GSEA', mustWork = TRUE),
     gs.db = "GSEA/GSEA_db/c2.cp.kegg.v7.0.symbols.gmt",
     
     
     output.directory = "GSEA/",
     doc.string= "kegg",
     #   non.interactive.run   = F,               
     reshuffling.type      = "sample.labels", 
     nperm                 = 1000,
     weighted.score.type   =  1,            
     nom.p.val.threshold   = -1,            
     fwer.p.val.threshold  = -1,            
     fdr.q.val.threshold   = 0.25,          
     topgs                 = 20,            
     adjust.FDR.q.val      = F,             
     gs.size.threshold.min = 1,             
     gs.size.threshold.max = 500,           
     reverse.sign          = F,             
     preproc.type          = 0,             
     random.seed           = 3338,          
     perm.type             = 0,             
     fraction              = 1.0,           
     replace               = F,             
     save.intermediate.results = F,         
     #   OLD.GSEA              = F,             
     use.fast.enrichment.routine = T, 
     
     collapse.dataset = TRUE, collapse.mode = 'max')


