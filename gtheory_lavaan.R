 ################################################
 ##   Univariate G-theory analysis using SEM   ## 
 ################################################

 ## Article: A Robust Indicator Mean-based Method 
  #         for Estimating Generalizability Theory Absolute Error Indices 
  #         within Structural Equation Modeling Frameworks



## Read me: 
 # Steps for analyzing one-, two- and three-facet generalizability theory (G-theory) designs: 
 # We have created code in R to simplify the process for analyzing the designs above.
 # Running these analyses primarily entails using the 'gtheory_lavaan' function we have created. 
 # Within the function call, you will identify the data source and specify details for the analysis. 
 # The data source must be arranged using the “wide format” in which a separate row is 
 # included for each person (or object of measurement), and a separate column is included for
 # each variable. Analysis details include specification of: (a) names of each variable analyzed, 
 # (b) names for each facet, (c) number of conditions for each facet, (d) estimation procedure, 
 # and (e) method of parameterization when using WLSMV estimation. Results are derived 
 # primarily from the lavaan package in R. Monte Carlo-based confidence intervals for G, global 
 # D coefficients, and variance components are obtained using the 'monteCarloCI' function from 
 # the semTools package. More specific details for running the analyses and an example follow.


## Function: gtheory_lavaan(data, facet, d_n, method, estimator, parameterization)

## Arguments:
 # data = wide format data containing the observed variables in the model. 
 # facet = Variable names indicating each facet. 
 #         Example: For a 2-facet design with variable names like "item1rater1", "item2rater1"..., 
 #         facet = c(f1 = "item", f2 = "rater")
 # d_n = number of conditions for each facet for the D-study. 
 # method = absolute error estimation method. Either "indicator_mean" or "Jorgensen"
 # estimator = SEM model estimation method such as "ULS", "WLSMV". See the lavaan function for alternative estimators.
 # parameterization = metric setting when "WLSMV" estimator is used. Either "delta" or "theta"
 
## Example:
 # G-theory analysis results (lavaan object)
  Results <- gtheory_lavaan(data= mydata, facet= c(f1="I", f2="O", f3="S"), 
                            d_n= c(12, 1, 4), method="indicator_mean", 
                            estimator="WLSMV", parameterization="theta")

 # Confidence intervals for G, D coefficients, and variance components
  library(semTools) 
  set.seed(123)
  monteCarloCI(Results, level=0.95) # level = numeric confidence level between 0-1
  


## gtheory_lavaan function code
library(lavaan) 

gtheory_lavaan<- function(data, facet, d_n, method="indicator_mean", estimator="ULS", parameterization="theta"){
  if(length(facet)==1) {
    colnames(data) <- gsub(paste0(facet["f1"], "(\\d+)"), "I\\1", sprintf("I%02d", 
                                                                          as.numeric(gsub(paste0(".*", facet["f1"], "(\\d+).*"), "\\1", colnames(data)))))
    name =colnames(data)
    
    ni_g =ncol(data)
    ni_d = d_n[1]
    ni = paste0("I",sprintf("%02d",1:ni_g))
    
    ### universe score and relative error
    model <-c(paste("person=~", paste0("1*", name, collapse = " + ")),
              paste0(name, "~~pi_var*", name),"person~~p_var*person",
              paste0(name,"~m_",name,"*1"),
              paste0("m:= (", paste0("m_", name, collapse = " + "), ")/", ni_g))
    
    cal <-c(
      paste0("Gcoef:= p_var/(p_var+pi_var/",ni_d,")"),
      paste0("Dcoef:= p_var/(p_var+pi_var/",ni_d,"+S_i/",ni_d,")"),
      "VC_p:= p_var",
      "VC_pXf1:= pi_var",
      "VC_f1:= S_i"
    )
    
    # WLSMV: threshold
    n_thre <-max(data)-min(data)
    thre <-c(paste(colnames(data),'|', paste(sapply(1:n_thre, function(i) paste0("Thr_", i, "*t", i)), collapse = "+")))
    
    resid1 <-'pi_var==1'
    
    if (method=="indicator_mean"){
      abs <-paste0("S_i:= (", paste0("(m_", name,"-m)^2", collapse = " + "), ")/(", ni_g, "-1)") 
      if (estimator=="ULS"){
        lavaan(model= c(model, abs, cal), orthogonal = T, data = data, estimator= "ULS")
        
      }else if(estimator=="WLSMV" & parameterization=="theta"){
        lavaan(model = c(model, resid1, abs, thre, cal), orthogonal = T, 
               data = data, ordered=name, estimator= "WLSMV", parameterization="theta")
      }
      else if(estimator=="WLSMV" & parameterization=="delta"){
        lavaan(model = c(model, abs, thre, cal), orthogonal = T, 
               data = data, ordered=name, estimator= "WLSMV", parameterization="delta")
      }
    } else if (method== "Jorgensen"){
      abs<-c("person~Mu*1",
             paste0("m_", name[ni_g], "==-1*(", paste0("m_", name, collapse = " + "), ")"),
             paste0("S_i:= (", paste0("(m_", name, "-m)^2", collapse = " + "), ")/(", ni_g, "-1)")
      )
      if (estimator=="ULS"){
        lavaan(model= c(model, abs, cal), orthogonal = T, data = data, estimator= "ULS")
        
      }else if(estimator=="WLSMV"){
        lavaan(model = c(model, resid1, abs, thre, cal), orthogonal = T, 
               data = data, ordered=name, estimator= "WLSMV", parameterization="theta")
      }
      
    }} else if (length(facet)==2) {
      colnames(data) <-gsub(paste0(facet["f1"], "(\\d+)",facet["f2"],"(\\d+)"), "O\\1I\\2",
                            sprintf("O%02dI%02d", 
                                    as.numeric(gsub(paste0(".*",facet["f2"], "(\\d+).*"), "\\1", colnames(data))),
                                    as.numeric(gsub(paste0(".*",facet["f1"], "(\\d+).*"), "\\1", colnames(data)))))
      name =colnames(data)
      
      ni_g =length(unique(substr(name, regexpr("I", name),regexpr("I", name)+2)))
      no_g =length(unique(substr(name, regexpr("O", name),regexpr("O", name)+2)))
      ni_d = d_n[1]
      no_d = d_n[2]
      
      ni = matrix(paste0("I",sprintf("%02d",1:ni_g))) 
      no = matrix(paste0("O",sprintf("%02d",1:no_g)))
      
      fact  <- function(combi) {
        if(ncol(combi)==2) {
          match1 <- grep(combi$f1, name, value=T)
          match2 <- grep(combi$f2, match1, value=T)
          return(paste0("p", combi$f1, combi$f2, "=~ ",
                        paste0("1*",match2, collapse = ' + ')))
        } else {
          return(paste0("p", combi, "=~ ", paste0("1*", grep(combi, name, value=T), collapse = ' + ')))
        }
      }
      
      model<-c(
        paste0("person=~", paste0("1*",name, collapse = ' + ')),
        sapply(1:nrow(ni), function(i) fact(ni[i,,drop=F])),
        sapply(1:nrow(no), function(i) fact(no[i,,drop=F])),
        "person~~p_var*person",
        paste0("p",ni,"~~pi_var*","p",ni),
        paste0("p",no,"~~po_var*","p",no),
        paste0(name,"~~pio_var*",name),
        paste0(name,"~m_",name,"*1"),
        paste0("m:=(", paste0("m_", name, collapse = ' +'),")/",
               ncol(data)))
      
      cal <-c(
        paste0("Gcoef:=p_var/(p_var+pi_var/",ni_d,"+po_var/",no_d,"+pio_var/",(ni_d*no_d),")"),
        paste0("Dcoef:= p_var/(p_var+pi_var/",ni_d,"+po_var/",no_d,"+pio_var/",(ni_d*no_d),
               "+S_i/",ni_d,"+S_o/",no_d,"+S_io/",(ni_d*no_d),")"),
        "VC_p:= p_var",
        "VC_pXf1:= pi_var",
        "VC_pXf2:= po_var",
        "VC_pXf1Xf2:= pio_var",
        "VC_f1:= S_i",
        "VC_f2:= S_o",
        "VC_f1Xf2:= S_io")
      
      # WLMSV
      n_thre <-max(data)-min(data)
      thre <-c(paste(colnames(data),'|', paste(sapply(1:n_thre, function(i) paste0("Thr_", i, "*t", i)), collapse = "+")))
      resid1 <-'pio_var==1'
      
      if (method=="indicator_mean"){
        meansum <- function(f, d) {
          if(length(f)==1){
            match <- grep(f, name, value = TRUE)
            paste0("m", f, ":=(", paste0("m_", match, collapse = 
                                           ' + '), ")/", d)
          } else { 
            match1 <- grep(f$f1, name, value=T)
            match2 <- grep(f$f2, match1, value=T)
            paste0("m", f$f1, f$f2, ":=(", paste0("m_", match2, collapse = ' + '), ")/", d)
          }
        }
        
        # absolute error: our procedure
        abs <-c(
          apply(ni, 1, function(i) meansum(i, no_g)),
          apply(no, 1, function(i) meansum(i, ni_g)),
          paste0("S_i:=(", paste0("(m", ni, "-m)^2", collapse = ' + '), ")/(", ni_g, "-1)"),
          paste0("S_o:=(", paste0("(m", no, "-m)^2", collapse = ' + '), ")/(", no_g, "-1)"),
          paste0("S_io:=(", paste0("(m_", name, "-m", substr(name, regexpr('I', name), regexpr('I', name) + 2),
                                   "-m", substr(name, regexpr('O', name), regexpr('O', name) + 2), "+m)^2", collapse = ' + '), ")/(", ni_g, "*", no_g, "-1)")
        )
        if (estimator=="ULS"){
          lavaan(model= c(model, abs, cal), orthogonal = T, data = data, estimator= "ULS")
        } else if (estimator=="WLSMV" & parameterization=="theta"){
          lavaan(model = c(model, resid1, abs, thre, cal), 
                 orthogonal = T, data = data, ordered=name, estimator=
                   "WLSMV", parameterization="theta")
        } else if (estimator=="WLSMV" & parameterization=="delta"){
          lavaan(model = c(model, abs, thre, cal), orthogonal = T, 
                 data = data, ordered=name, estimator= "WLSMV", parameterization="delta")
        }}
      else if (method=="Jorgensen"){
        sumzero <- function(f) {
          if (length(f)==1) {paste0("(", paste0("m_", grep(f, name, value = T), collapse =' + '), ")==0")
          } else {
            match1 <- grep(f$f1, name, value=T)
            match2 <- grep(f$f2, match1, value=T)
            paste0("(", paste0("m_", match2, collapse =' + '), ")==0")}}
        paste0("m_", name[ncol(data)], " == -1*(", paste0("m_", name, collapse= ' + '), ")")
        
        abs <-c("person~Mu*1", 
                paste0("p", ni, "~mu_", ni, "*1"),
                paste0("p", no, "~mu_", no, "*1"),
                paste0("m_", colnames(data)[ncol(data)], "== -1*(",
                       paste0("m_", colnames(data), collapse = ' + '), ")"),
                paste0("(", paste0("mu_", ni, collapse = ' + '), ")==0"),
                paste0("(", paste0("mu_", no, collapse = ' + '), ")==0"),
                apply(ni, 1, function(f) sumzero(f)),
                apply(no, 1, function(f) sumzero(f)),
                paste0("S_i:=(", paste0("(mu_", ni, ")^2", collapse = ' + '), ")/(", ni_g, "-1)"),
                paste0("S_o:=(", paste0("(mu_", no, ")^2", collapse = ' + '), ")/(", no_g,"-1)"),
                paste0("S_io:=(", paste0("(m_", colnames(data), ")^2", collapse = ' + '), ")/((", ni_g, "-1)*(", no_g,"-1))")
        )
        if (estimator=="ULS") {lavaan(model= c(model, abs, cal),
                                      orthogonal = T, data = data, estimator= "ULS")
        } else if (estimator=="WLSMV") {
          lavaan(model = c(model, resid1, abs, thre, cal), 
                 orthogonal = T, data = data, ordered=name, 
                 estimator= "WLSMV", parameterization="theta")
        }
      }
    } else if(length(facet)==3){
      colnames(data)<-gsub(paste0(facet["f3"],"(\\d+)",facet["f2"], "(\\d+)",
                                  facet["f1"], "(\\d+)"), "O\\1S\\2I\\3",
                           sprintf("O%02dS%02dI%02d", 
                                   as.numeric(gsub(paste0(".*",facet["f2"],  "(\\d+).*"), "\\1", colnames(data))), 
                                   as.numeric(gsub(paste0(".*",facet["f3"],  "(\\d+).*"), "\\1", colnames(data))), 
                                   as.numeric(gsub(paste0(".*",facet["f1"], "(\\d+).*"), "\\1", colnames(data)))))
      name = colnames(data)
      
      ni_g =length(unique(substr(name, regexpr("I", name),regexpr("I", name)+2)))
      no_g =length(unique(substr(name, regexpr("O", name),regexpr("O", name)+2)))
      ns_g =length(unique(substr(name, regexpr("S", name),regexpr("S", name)+2)))
      ni_d = d_n[1]
      no_d = d_n[2]
      ns_d = d_n[3]
      
      ni = matrix(paste0("I",sprintf("%02d",1:ni_g)))
      no = matrix(paste0("O",sprintf("%02d",1:no_g)))
      ns = matrix(paste0("S",sprintf("%02d",1:ns_g)))
      nio = as.data.frame(matrix(unlist(expand.grid(ni, no)), ncol=2, dimnames=list(NULL,c("f1","f2"))))
      nis = as.data.frame(matrix(unlist(expand.grid(ni, ns)), ncol=2, dimnames=list(NULL,c("f1","f2"))))
      nos = as.data.frame(matrix(unlist(expand.grid(no, ns)), ncol=2, dimnames=list(NULL,c("f1","f2"))))
      
      fact  <- function(combi) {
        if(ncol(combi)==2) {
          match1 <- grep(combi$f1, name, value=T)
          match2 <- grep(combi$f2, match1, value=T)
          return(paste0("p", combi$f1, combi$f2, "=~ ",
                        paste0("1*",match2, collapse = ' + ')))
        } else {
          return(paste0("p", combi, "=~ ", paste0("1*", grep(combi, name, value=T), collapse = ' + ')))
        }
      }
      
      model <- c(
        paste0("person=~", paste0("1*",name, collapse = ' + ')),
        sapply(1:nrow(ni), function(i) fact(ni[i,,drop=F])),    
        sapply(1:nrow(no), function(i) fact(no[i,,drop=F])),   
        sapply(1:nrow(ns), function(i) fact(ns[i,,drop=F])),   
        sapply(1:nrow(nio), function(i) fact(nio[i, ])),
        sapply(1:nrow(nis), function(i) fact(nis[i, ])),
        sapply(1:nrow(nos), function(i) fact(nos[i, ])),
        "person~~p_var*person",
        paste0("p", ni, "~~pi_var*", "p", ni),
        paste0("p", no, "~~po_var*", "p", no),
        paste0("p", ns, "~~ps_var*", "p", ns),
        paste0("p", paste0(nio$f1,nio$f2), 
               "~~pio_var*", "p", paste0(nio$f1,nio$f2)),
        paste0("p", paste0(nis$f1,nis$f2), "~~pis_var*", "p",
               paste0(nis$f1,nis$f2)),
        paste0("p", paste0(nos$f1,nos$f2), "~~pos_var*", "p",
               paste0(nos$f1,nos$f2)),
        paste0(name, "~~pios_var*", name),
        paste0(name, "~m_", name, "*1"),
        paste0("m:=(", paste("m_", name, collapse = ' + '), ")/", ncol(data)))
      
      cal <-c(
        paste0("Gcoef:= p_var/(p_var+pi_var/",ni_d,"+po_var/",no_d,"+ps_var/",ns_d,
               "+pio_var/",(ni_d*no_d),"+pis_var/",(ni_d*ns_d),"+pos_var/",
               (no_d*ns_d),"+pios_var/",(ni_d*ns_d*no_d),")"),
        paste0("Dcoef:= p_var/(p_var+pi_var/",ni_d,"+po_var/",no_d,"+ps_var/",ns_d,
               "+pio_var/",(ni_d*no_d),"+pis_var/",(ni_d*ns_d),"+pos_var/",(no_d*ns_d),
               "+pios_var/",(ni_d*ns_d*no_d),"+S_i/",ni_d,"+S_o/",no_d,"+S_s/",ns_d,
               "+S_io/",(ni_d*no_d),"+S_is/",(ni_d*ns_d),"+S_os/",(no_d*ns_d),
               "+S_ios/",(ni_d*no_d*ns_d),")"),
        "VC_p:= p_var",
        "VC_pXf1:= pi_var",
        "VC_pXf2:= po_var",
        "VC_pXf3:= ps_var",
        "VC_pXf1Xf2:= pio_var",
        "VC_pXf1Xf3:= pis_var",
        "VC_pXf2Xf3:= pos_var",
        "VC_pXf1Xf2Xf3:= pios_var",
        "VC_f1:= S_i",
        "VC_f2:= S_o",
        "VC_f3:= S_s",
        "VC_f1Xf2:= S_io",
        "VC_f1Xf3:= S_is",
        "VC_f2Xf3:= S_os",
        "VC_f1Xf2Xf3:= S_ios")
      
      # WLMSV
      resid1 <-'pios_var==1'
      n_thre <-max(data)-min(data)
      thre <-c(paste(colnames(data),'|', paste(sapply(1:n_thre, function(i)
        paste0("Thr_", i, "*t", i)), collapse = "+")))
      
      if (method=="indicator_mean"){
        meansum <- function(f, d) {
          if(length(f)==1){
            match <- grep(f, name, value = TRUE)
            paste0("m", f, ":=(", paste0("m_", match, collapse = 
                                           ' + '), ")/", d)
          } else { 
            match1 <- grep(f$f1, name, value=T)
            match2 <- grep(f$f2, match1, value=T)
            paste0("m", f$f1, f$f2, ":=(", paste0("m_", match2, collapse = ' + '), ")/", d)
          }
        }
        
        name_s<-function(f, name){
          substr(name,regexpr(f,name),regexpr(f,name)+2)}
        
        abs <-c(
          apply(ni, 1, function(i) meansum(i, no_g*ns_g)),
          apply(no, 1, function(i) meansum(i, ni_g*ns_g)),
          apply(ns, 1, function(i) meansum(i, ni_g*no_g)),
          apply(nio, 1, function(r){meansum(setNames(
            as.list(r), c("f1", "f2")), ns_g)}),
          apply(nis, 1, function(r){meansum(setNames(
            as.list(r), c("f1", "f2")), no_g)}),
          apply(nos, 1, function(r) {meansum(setNames(
            as.list(r), c("f1", "f2")), ni_g)}),
          paste0("S_i:=(",paste0("(m", ni, "-m)^2", 
                                 collapse = ' + '), ")/", ni_g-1),
          paste0("S_o:=(",paste0("(m", no, "-m)^2", collapse = 
                                   ' + '), ")/", no_g-1),
          paste0("S_s:=(",paste0("(m", ns, "-m)^2", collapse = 
                                   ' + '), ")/", ns_g-1),
          paste0("S_io:=(",paste0("(m", paste0(nio$f1,nio$f2),
                                  " - m",nio$f1," - m",nio$f2, " + m)^2", 
                                  collapse = ' + '), ")/", ni_g*no_g-1),
          paste0("S_is:=(",paste0("(m", paste0(nis$f1,nis$f2)," - m",nis$f1,
                                  " - m",nis$f2, " + m)^2", collapse = ' + '),
                 ")/", ni_g*ns_g-1),
          paste0("S_os:=(",paste0("(m", paste0(nos$f1,nos$f2)," - m",nos$f1,
                                  " - m",nos$f2, " + m)^2", collapse = ' + '), 
                 ")/", no_g*ns_g-1),
          paste0("S_ios:=(",capture.output(cat(paste0("(m_",name,"-m",
                                                      name_s("I",name),name_s("O",name),"-m",
                                                      name_s("I",name),name_s("S",name),"-m",
                                                      name_s("O",name),name_s("S",name),"+m",
                                                      name_s("I",name),"+m",name_s("O",name),"+m",
                                                      name_s("S",name),"-m)^2"), sep=' + ')),")/",
                 (ni_g*ns_g*no_g-1), collapse=""))
        
        if (estimator=="ULS"){
          lavaan(model= c(model, abs, cal), orthogonal = T, data = data, estimator= "ULS")
        } else if (estimator=="WLSMV" & parameterization=="theta"){
          lavaan(model = c(model, resid1, abs, thre, cal), orthogonal = T, data = data,
                 ordered=name, estimator= "WLSMV", parameterization="theta")
        } else if (estimator=="WLSMV" & parameterization=="delta"){
          lavaan(model = c(model, abs, thre, cal), orthogonal = T, data = data,
                 ordered=name, estimator= "WLSMV", parameterization="delta")
        }
      } else if (method=='Jorgensen'){
        abs<-c("person~Mu*1", 
               paste0("p",ni,"~m_",ni,"*1"),
               paste0("p",no,"~m_",no,"*1"),
               paste0("p",ns,"~m_",ns,"*1"),
               paste0("p",rep(ni,each=no_g),no,"~m_",rep(ni,each=no_g),no,"*1"),
               paste0("p",rep(ni,each=ns_g),ns,"~m_",rep(ni,each=ns_g),ns,"*1"),
               paste0("p",rep(no,each=ns_g),ns,"~m_",rep(no,each=ns_g),ns,"*1"),
               paste0("m_", name[ncol(data)], " == -1*(", paste0("m_", name, collapse= ' + '), ")"),
               apply(ni, 1, function(f) sumzero(f)),
               apply(no, 1, function(f) sumzero(f)),
               apply(ns, 1, function(f) sumzero(f)),
               apply(nio, 1, function(f) {sumzero(setNames(as.list(f), c("f1", "f2")))}),
               apply(nis, 1, function(f) {sumzero(setNames(as.list(f), c("f1", "f2")))}),
               apply(nos, 1, function(f) {sumzero(setNames(as.list(f), c("f1", "f2")))}),
               paste0("(",paste0("m_",ni,collapse="+"),")==0"),
               paste0("(",paste0("m_",no,collapse="+"),")==0"),
               paste0("(",paste0("m_",ns,collapse="+"),")==0"),
               paste0("(",paste0("m_",paste0(nio$f1,nio$f2),collapse="+"),")==0"),
               paste0("(",paste0("m_",paste0(nis$f1,nis$f2),collapse="+"),")==0"),
               paste0("(",paste0("m_",paste0(nos$f1,nos$f2),collapse="+"),")==0"),
               paste0("S_i:=(",paste0("(m_", ni,")^2", collapse = ' + '), ")/", ni_g-1),
               paste0("S_o:=(",paste0("(m_", no,")^2", collapse = ' + '), ")/", no_g-1),
               paste0("S_s:=(",paste0("(m_", ns,")^2", collapse = ' + '), ")/", ns_g-1),
               paste0("S_io:=(",paste0("(m_", paste0(nio$f1,nio$f2),")^2", collapse = ' + '),
                      ")/", ni_g*no_g-1),
               paste0("S_is:=(",paste0("(m_", paste0(nis$f1,nis$f2),")^2", collapse = ' + '),
                      ")/", ni_g*ns_g-1),
               paste0("S_os:=(",paste0("(m_", paste0(nos$f1,nos$f2),")^2", collapse = ' + '),
                      ")/", no_g*ns_g-1),
               paste0("S_ios:=(",paste0("(m_", name,")^2", collapse = ' + '), ")/", ni_g*no_g*ns_g-1)
        )
        
        if (estimator=="ULS"){
          lavaan(model= c(model, abs, cal), orthogonal = T, data = data, estimator= "ULS")
        } else if (estimator=="WLSMV"){
          lavaan(model = c(model, resid1, abs, thre, cal), orthogonal = T, 
                 data = data, ordered=name, estimator= "WLSMV", parameterization="theta")
        }
      }
    }
}
