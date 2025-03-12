#Gérer 'ehahelper' (recheker avec Fred)

# Load packages
library(lubridate)
library(tidyverse)
library(nlme)
library(lme4)
library(ggplot2)
library(coxme)
library(boot)
#library(ehahelper)
library(ggridges)
library(ggnewscale)
library(sf)
library(raster)
library(mapview)

# —————————————————————————————————————————————————————————————————————————————————————
# Analysis 1. Comparison of time spent foraging by cort- and placebo-treated geese ####
# —————————————————————————————————————————————————————————————————————————————————————

## Data formating
MSA_raw <- readRDS("behavior_data_IAO_2021.RDS")

## Regrouping clusters into behaviors:
MSA_behavior <- MSA_raw%>%
  mutate(foraging=case_when(
    EMcluster5==3~1,
    EMcluster5==2~0,
    TRUE~0),
    forandother= case_when(
      EMcluster5==3~1,
      EMcluster5==2~1,
      TRUE~0),
    other= case_when(
      EMcluster5==3~0,
      EMcluster5==2~1,
      TRUE~0))

#Removing pair number 9 because HJ was injured during the capture
MSA_day<-MSA_behavior%>%
  filter(paire!=9)

#Add number of observation per ID, day, periode
data_MSA<- MSA_day%>%
  group_by(ID,JT,periode, trt, paire,numero_capture)%>%
  summarise(nb_obs= n(),
            prop_foraging = sum(foraging)/n(),
            prop_all= sum(forandother)/n(),
            prop_other = sum(other)/n())


#Filter, relevel and clean dataset
data_MSA <- data_MSA[,c("prop_foraging","prop_all", "prop_other", "trt", "periode", "JT", "paire","numero_capture", "ID", "nb_obs")]
data_MSA <- data_MSA %>% mutate(trt=factor(trt, levels = c("placebo", "cort")),
                                periode = factor(periode, levels = c("day", "evening", "morning")),
                                paire = factor(paire),
                                ID = factor(ID))

#................................................................................................
### Table 1. Generalized linear mixed model (GLMM) for the proportion of time spent foraging ####
#................................................................................................

## Model fitting
GLM1 <- MASS::glmmPQL(prop_foraging~ trt + JT + periode + trt:JT, random = list(~1|paire, ~1|ID), family = quasibinomial, weight=nb_obs, data=data_MSA[data_MSA$JT>1 ,], na.action = na.omit)
summary(GLM1)

#Get confidence intervals for all coefficient estimates
GLM1.int<-as_tibble(intervals(GLM1)$fixed)%>%
  mutate(par=rownames(intervals(GLM1)$fixed))
GLM1.int

#...................................
### Predictions and effect size ####
#...................................

# get model predictions

# build newdata
newdata <- data.frame(expand.grid(trt = levels(data_MSA$trt),
                                  periode= levels(data_MSA$periode),
                                  JT = seq(from= min(data_MSA$JT), to= max(data_MSA$JT), by=1)))
                      
newdata$trt <- factor(newdata$trt, levels=c("placebo", "cort")) 

# model matrix for newdata
X_new <- model.matrix(~trt + JT + periode + trt:JT, data = newdata)

# fixed effect estimates
beta <- fixef(GLM1)

# logit predictions
newdata$pred_logit <- X_new %*% beta

# standard errors
vcov_mat <- vcov(GLM1) # Extract variance-covariance matrix
se_logit <- sqrt(diag(X_new %*% vcov_mat %*% t(X_new))) 

# 95% confidence intervals on logit scale
newdata$logit_low <- newdata$pred_logit - 1.96 * se_logit
newdata$logit_high <- newdata$pred_logit + 1.96 * se_logit

# probability scale
newdata$pred_prob <- plogis(newdata$pred_logit)
newdata$prob_low <- plogis(newdata$logit_low)
newdata$prob_high <- plogis(newdata$logit_high)

# filter to keep only day for graph
newdata_plot <- newdata[newdata$periode=='day' & newdata$JT>=2,] 

## Effect size foraging rate
desired_JT=2:10
forage_cort <- newdata_plot[newdata_plot$trt=='cort'&newdata_plot$JT%in%desired_JT,"pred_prob"][,1]
forage_placebo <- newdata_plot[newdata_plot$trt=='placebo'&newdata_plot$JT%in%desired_JT,"pred_prob"][,1]

data.frame(JT=desired_JT, 
           absolute_change=forage_cort-forage_placebo,
           relative_change=(forage_cort-forage_placebo)/forage_cort)

# overall predicted increase in foraging time
auc_cort <- sum(forage_cort * diff(c(1, desired_JT)))
auc_placebo <- sum(forage_placebo  * diff(c(1, desired_JT)))
auc_diff <- auc_cort - auc_placebo
percentage_increase <- (auc_diff / auc_placebo) * 100
percentage_increase


#.......................................................
### Fig.1. Foraging intensity of greater snow goose ####
#.......................................................

library(ggridges)
library(ggnewscale)

data_MSA_sum <- data_MSA %>% 
  filter(periode=="day") %>% 
  group_by(JT, trt) %>% 
  summarize(prop_forage = mean(prop_foraging),
            se= sd(prop_foraging)/sqrt(n())) %>% 
  mutate(prop_low=prop_forage-1.96*se,
         prop_high=prop_forage+1.96*se
  )

dodge_width <- 0.2
newdata_plot$trt <- factor(newdata_plot$trt, levels=c("cort", "placebo")) 
data_MSA_day <- data_MSA[data_MSA$periode == "day",]
data_MSA_day$trt <- factor(data_MSA_day$trt, levels=c("cort", "placebo"))
data_MSA_day$group <- ifelse(data_MSA_day$trt=="cort", paste0("B",data_MSA_day$JT), paste0("A",data_MSA_day$JT))

A <- ggplot(data_MSA_day, aes(y= JT, x= prop_foraging, fill = trt)) +
  
  # line conf intervals
  geom_line(data=newdata_plot, aes(y=JT-dodge_width, x=prob_low, col=trt), show.legend=F, lty=2, alpha=1)+
  geom_line(data=newdata_plot, aes(y=JT-dodge_width, x=prob_high, col=trt), show.legend=F, lty=2, alpha=1)+
  
  # points and errorbar
  geom_errorbar(data=data_MSA_sum, aes(y=JT-dodge_width, x=prop_forage, xmin= prop_low, xmax=prop_high),
                show.legend = F, alpha=1, width=0, size=1)+
  
  geom_point(data = subset(data_MSA_sum, trt == "placebo"), aes(y = JT - dodge_width, x = prop_forage, fill = trt),
             show.legend = F, alpha = 1, size = 4.5, stroke = 1.7, pch = 21) +
  
  geom_point(data = subset(data_MSA_sum, trt == "cort"), aes(y = JT - dodge_width, x = prop_forage, fill = trt),
             show.legend = F, alpha = 1, size = 4.5, stroke = 1.7, pch = 21) +
  
  # lines predictions
  geom_line(data=newdata_plot, aes(y=JT-dodge_width, x=pred_prob, col=trt), show.legend=F, lty=1, lwd=1.2, alpha=1)+
  
  scale_fill_manual(values = c("#A9CDDB", "#F3977F"), breaks = c("placebo", "cort"), labels= c("PLACEBO", "CORT"), guide = "none")+
  
  # new scale to display violin plot
  new_scale_fill() +
  geom_density_ridges2(data=data_MSA_day, aes(y= JT, x= prop_foraging, fill= trt, group=group),
                       scale =0.7, show.legend = T, alpha=0.4, bandwidth = 0.045, col= NA)+
  
  scale_y_continuous(n.breaks=10)+
  coord_flip(xlim= c(0,0.8))+
  xlab("Proportion of time foraging")+
  ylab("Days after pellet implantation")+
  theme_minimal()+ 
  scale_color_manual(values = c("#F3977F", "#A9CDDB"), breaks = c("cort", "placebo"), labels= c("CORT", "PLACEBO"))+
  scale_fill_manual(values=c("#E73000", "#529CB5"), breaks = c("cort", "placebo"), labels= c("CORT", "PLACEBO"),
                    guide = guide_legend(override.aes = list(color = NA)))+
  theme(panel.grid=element_blank(),
        legend.title=element_blank(),
        legend.text = element_text(size=14, colour="black"),
        axis.text=element_text(size=17, colour="black"),
        axis.title=element_text(margin =margin(0,10,20,), size=18,face="bold"),
        axis.line = element_line(color = 'black'),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 8)))
A

ggsave("Figure_foraging.png", A, width = 10, height = 5.8, dpi=300, bg = 'white')
ggsave("Figure_1.pdf", A, width = 10, height = 5.8, dpi=300, bg = 'white')

#....................................................
### R2 computation for time spent foraging model ####
#....................................................

#Residual variance and quantities to compute R2 statistic are not straightforward with quasibinomial distributions.

#The square of the dispersion parameter (residual SD) can be thought of as an approximation of residual variance in the model. Thus, if 
#this quantity is substracted from that of the same model without the fixed effects it should approximate the variance explained by the 
#fixed effects in GLM1 (and be used to approximate a marginal R2). Similarly, the difference between the dispersion parameter of GLM1 
#and a that of a model without any fixed or random effects should approximate the variance explained by the model's fixed AND random 
#effects (and be used to approximate a conditional R2). The proportion of variance "explained" by the random and fixed effects of the 
#model can be approximated using this quantity using the method of Nagakawa and Schielzeth (2013).

#Create dummy variable as grouping factor to use in glmmPQL for model without random effects
data_MSA2<-data_MSA
data_MSA2$dummy<-0

#Model without fixed effects
GLM0 <- MASS::glmmPQL(prop_foraging~ 1, random = list(~1|paire, ~1|ID), family = quasibinomial, weight=nb_obs, data=data_MSA[data_MSA$JT>1 ,], na.action = na.omit)

#Variance explained by 'paire' near 0 and very similar for ID in both models. OK to approximate marginal R2
VarCorr(GLM1)
VarCorr(GLM0)

#Model without fixed or random effects
GLM00 <- MASS::glmmPQL(prop_foraging~ 1, random = list(~1|dummy), family = quasibinomial, weight=nb_obs, data=data_MSA2[data_MSA2$JT>1 ,], na.action = na.omit)
VarCorr(GLM00)[1,'StdDev'] #Virtually no variance explained by the 'dummy' random effect, OK to approximate conditional R2

#Residual variance of each model
resid_var_m1 <-summary(GLM1 )$sigma^2-1 #residual variance of model
resid_var_m0 <-summary(GLM0 )$sigma^2-1 #residual variance when no fixed effects fitted
resid_var_m00<-summary(GLM00)$sigma^2-1 #residual variance when no fixed or random effects fitted

#Total variance is residual variance of model without any fixed or random effects
tot_var <- resid_var_m00
#Variance explained by random effects: total variance - residual variance random effects only
var_ranef <- tot_var - resid_var_m0
#Variance explained by the fixed effects: 
#residual variance no effects - residual variance full model - variance explained by random effects
var_fixed <- tot_var - resid_var_m1 - var_ranef

#Marginal R2
(R2_marginal <- var_fixed/tot_var)
#Conditional R2
(R2_conditional <- (var_fixed + var_ranef)/tot_var)

#....................................................................................
### Appendix S1: Table S1 - Analysis with foraging and mixed activities combined ####
#....................................................................................

GLMS1 <- MASS::glmmPQL(prop_all~ periode + trt*JT, random = list(~1|paire, ~1|ID), family = quasibinomial, weight=nb_obs, data=data_MSA[data_MSA$JT>1 ,], na.action = na.omit)
summary(GLMS1)
intervals(GLMS1)

#Same steps to compute R2
GLMS0 <- MASS::glmmPQL(prop_all~ 1, random = list(~1|paire, ~1|ID), family = quasibinomial, weight=nb_obs, data=data_MSA[data_MSA$JT>1 ,], na.action = na.omit)

#In both models, variance explained by 'paire' near 0 and very similar for ID. OK to approximate marginal R2
VarCorr(GLMS1)
VarCorr(GLMS0)

#Model without fixed or random effects
GLMS00 <- MASS::glmmPQL(prop_all~ 1, random = list(~1|dummy), family = quasibinomial, weight=nb_obs, data=data_MSA2[data_MSA2$JT>1 ,], na.action = na.omit)
VarCorr(GLMS00)[1,'StdDev'] #No variance explained by the 'dummy' random effect, OK to approximate conditional R2

#Residual variance of each model
resid_var_m1 <-summary(GLMS1 )$sigma^2-1 #residual variance of model
resid_var_m0 <-summary(GLMS0 )$sigma^2-1 #residual variance when no fixed effects fitted
resid_var_m00<-summary(GLMS00)$sigma^2-1 #residual variance when no fixed or random effects fitted

#Total variance is residual variance of model without any fixed or random effects
tot_var <- resid_var_m00
#Variance explained by random effects: total variance - residual variance random effects only
var_ranef <- tot_var - resid_var_m0
#Variance explained by the fixed effects: 
#residual variance no effects - residual variance full model - variance explained by random effects
var_fixed <- tot_var - resid_var_m1 - var_ranef

#Marginal R2
(R2_marginal <- var_fixed/tot_var)
#Conditional R2
(R2_conditional <- (var_fixed + var_ranef)/tot_var)



# —————————————————————————————————————————————————————————————————————————————————————
# Analysis 2. Comparison of departure dates between cort- and placebo-treated geese ####
# —————————————————————————————————————————————————————————————————————————————————————

# load data
a <- readRDS("departure_dates.RDS")
# Filter pair for which one bird either died or had a defective collar (we never consider pair no. 9 because bird was injured during capture)
'%!in%' <- function(x,y)!('%in%'(x,y))
a<- a[a$paire %!in% c(9,10,11,22,25,26), ]

# Format dataframe
a$depart_date<- as.numeric(format(as.POSIXct(strptime(a$depart_date , "%Y-%m-%d")),"%d"))
a$departure_date<- as.numeric(a$julian_departure)
a$date_deploiement<- as.numeric(format(as.POSIXct(strptime(a$Deploy.On.Timestamp , "%Y-%m-%d")),"%j"))
a$trt<- factor(a$trt, levels = c("cort", "placebo"))
a$date_capt<- a$date_deploiement-(min(a$date_deploiement)-1)

#.....................................
### Compute median departure date ####
#.....................................

# Define a bootstrapping function for the median
bootstrap_median <- function(data, indices) {
  sample_data <- data[indices]
  return(median(sample_data))
}

# Apply bootstrapping to calculate median and CI for each group
#install.packages("boot")
library(boot)
result <- a %>%
  group_by(trt) %>%
  do({
    boot_result <- boot(
      data = .$depart_date, 
      statistic = bootstrap_median, 
      R = 1000  # Number of bootstrap samples
    )
    boot_ci <- boot.ci(boot_result, type = "perc")  # Percentile CI
    
    data.frame(
      trt = unique(.$trt),
      median = median(.$depart_date),
      CI_L = boot_ci$percent[4],  # 2.5th percentile
      CI_U = boot_ci$percent[5]   # 97.5th percentile
    )
  })
print(result)

#..................................................
### Table 2. Departure date model (cox hazard) ####
#..................................................

# Format dataframe to one line per day with a column is_departed (0/1)
instance_df <- data.frame()
for(i in c(1:nrow(a))){
  current_line <- a[i,]
  time_tracked <- current_line$departure_date-current_line$date_deploiement
  
  all_instance_lines <- current_line[rep(1, each = time_tracked), ]
  all_instance_lines$day_since_trt = c(1:time_tracked)
  all_instance_lines$current_day = all_instance_lines$date_deploiement + all_instance_lines$day_since_trt
  
  all_instance_lines$is_departed <- with(all_instance_lines, case_when(current_day < departure_date ~ 0,
                                                                       current_day == departure_date ~ 1,
                                                                       current_day > departure_date  ~ NA))
  instance_df <- rbind(instance_df, all_instance_lines)
}

# Fit cox model
library(coxme)
cme <- coxme(Surv(day_since_trt, is_departed) ~ trt*date_capt + (1|paire), data= instance_df, x=T)
summary(cme)

#Get table with confidence intervals
(est_plc <- as.numeric(exp(cme$coefficients)))
(CI <- confint(cme))
coef_table<-as.data.frame(cbind(cme$coefficients, CI))
colnames(coef_table)<-c('Estimate','LowerCI','UpperCI')
coef_table

#...................................................................................
### Figure 2A. Distribution of departure date for cort and placebo individuals #####
#...................................................................................

library(ggridges)
a<- a%>%mutate(trt=factor(trt, levels = c("placebo", "cort")))

A<- ggplot()+
  ggridges::geom_density_ridges2(mapping= aes(y=trt, x=depart_date,fill= trt), data= a, scale=0.9, bandwidth =0.9)+
  geom_dotplot(aes(y=trt, x=depart_date,fill= trt), data= a[a$trt=="placebo",], position = position_nudge(y=0), fill="gray45")+
  geom_dotplot(aes(y=trt, x=depart_date,fill= trt), data= a[a$trt=="cort",], position = position_nudge(y=0), fill="gray45")+
  
  scale_x_continuous(n.breaks = 10)+
  scale_y_discrete(labels = c("PLACEBO", "CORT"))+
  scale_fill_manual(values = c("#F3977F", "#A9CDDB"), breaks = c("cort", "placebo"))+
  
  geom_errorbar(mapping = aes(y= trt, xmin= CI_L, xmax= CI_U), data = result, width= 0.12, size=0.7, position= position_nudge(x = 0, y = -0.11))+
  geom_point(mapping = aes(y=trt, x=median, fill= trt),data= result, size=4.6, position= position_nudge(x = 0, y = -0.11), pch=21, stroke=1.2)+
  
  
  labs(y= "", x = "Departure date in May")+
  
  theme_classic()+
  theme(text = element_text(size = 20),
        legend.position = "none")
A

#......................................................................................
### Figure 2B. Predicted departure probabilities for cort and placebo individuals #####
#......................................................................................

# for-loop to generate predictions from the cox model
selectedcapdates<-sort(unique(instance_df$date_capt))

#We use the predict_coxme function from the ehahelper package. This package must be installed from source
#using devtools. To simplify we simply reproduce the needed function here directly from 
#the package's github page: https://github.com/junkka/ehahelper/blob/master/R/predict_coxme.R
predict_coxme <- function(object, 
                          newdata = NULL, 
                          type = c("lp", "risk"), 
                          se.fit = FALSE,
                          strata_ref = TRUE){
  
  if (!inherits(object, 'coxme'))
    stop("Primary argument much be a coxme object")
  
  type <- match.arg(type)
  n <- object$n[2]
  Terms <- delete.response(terms(object))
  has_strata <- !is.null(attr(Terms, "specials")$strata) 
  if (has_strata) 
    has_strata <- ifelse(length(attr(Terms, "specials")$strata) == 0, FALSE, has_strata)
  has_newdata  <- !is.null(newdata)
  
  if (!se.fit & type == "lp" & !has_newdata) return(object$linear.predictor)
  
  coef <- fixed.effects(object)
  mf <- survival:::model.frame.coxph(object)
  
  # boot.ci
  
  
  
  if (has_newdata){
    m <- model.frame(Terms, newdata)
  } else {
    m <- mf
  }
  
  # if strata update terms
  if (has_strata){
    strata_terms <- untangle.specials(Terms, "strata")
    Terms2 <- Terms[-strata_terms$terms]
  } else {
    Terms2 <- Terms
  }
  
  if (has_newdata){
    mm <- model.matrix(Terms2, m)
    mm <- mm[ ,-1]
  }
  
  # has strata and reference is strata
  # calculate strata means
  if (has_strata & strata_ref){
    # Full model matrix
    x <- model.matrix(Terms, data = mf)
    
    oldstrat <- mf[[strata_terms$vars]]
    xmeans <- rowsum(x, oldstrat)/as.vector(table(oldstrat))
  }
  
  if (!has_newdata){
    # extract all cols in x which matches Terms
    mm <- model.matrix(Terms2, data =mf)[ ,-1]
    m <- mf
  }
  
  if (has_strata & strata_ref){
    newstrat <- m[[strata_terms$vars]]
    mm <- mm - xmeans[match(newstrat, row.names(xmeans)), colnames(mm)]
  } else {
    mm <- mm - rep(object$means, each = nrow(m))
  }
  
  # if (!has_newdata & !has_strata){
  #   pred <- object$linear.predictor
  # }
  if (length(coef) == 1){
    pred <- mm * coef
  } else {
    pred <- (mm %*% coef)
  }
  
  
  if (se.fit) se <- sqrt(rowSums((mm %*% vcov(object)) * mm))
  if (type == "risk"){
    pred <- exp(pred)
    if (se.fit) se <- se * sqrt(pred)
  }
  if (se.fit) list(fit = pred, se.fit = se)
  else pred
}

for(capdat in selectedcapdates){
  # build dummy dataset
  pred_data = data.frame(expand.grid(trt=unique(instance_df$trt),
                                     date_capt=capdat))
  #get predictions
  #library(ehahelper)
  preds <- predict_coxme(cme, pred_data, se.fit=TRUE, type = "lp")
  pred_data$est <- preds$fit
  pred_data$upperCI <- preds$fit + 1.96*preds$se.fit
  pred_data$lowerCI <- preds$fit - 1.96*preds$se.fit
  
  
  #build baseline hazard for mixed-effect cox model
  breslow_est_adj_inter <- function(time, status, X, B, fit_cox){
    data <- data.frame(time,status,X)
    data <- data[order(data$time), ]
    t    <- unique(data$time)
    k    <- length(t)
    h    <- rep(0,k)
   
    for(i in 1:k) {
      
      LP_sample <- sum(fit_cox$means * coef(fit_cox)) 
      
      #Individual linear predictor for interaction model
      LP_indiv <- (coef(fit_cox)['trtplacebo']*(data$trtplacebo))+(coef(fit_cox)['date_capt']*(data$date_capt))+(coef(fit_cox)['trtplacebo:date_capt']*(data$trtplacebo)*data$date_capt) 
      
      lp_centered <- (LP_indiv - LP_sample)[data$time>=t[i]]
      risk <- exp(lp_centered)
      h[i] <- sum(data$status[data$time==t[i]]) / sum(risk)
    }
    
    res <- cumsum(h)
    return(res)
  }
  
  cme <- coxme(Surv(day_since_trt, is_departed) ~ trt*date_capt + (1|paire), data= instance_df, x=T)
  H0_basehaz_cme <- cbind(breslow_cme=breslow_est_adj_inter(time=instance_df$day_since_trt, status=instance_df$is_departed, X=cme$x, B=cme$coefficients, fit_cox = cme), time=unique(instance_df$day_since_trt))
  
  
  # building dataset with confidence interval
  bh <- as_tibble(H0_basehaz_cme)
  colnames(bh) <- c("hazard","day_since_trt")
  
  #Trick to make CI appear until last point because geom_step() only plots a line until the last point
  bh<-rbind(bh,bh[22,])
  bh[nrow(bh),2]<-nrow(bh)
  bh_pred <- rbind(bh, bh)
  bh_pred$trt <- rep(c("placebo", "cort"), each = nrow(bh))
  
  ### calculate coxme estimate using the model function: h(t)=h0(t)exp(b1X1) with estimates obtained with coxme models
  bh_pred <- bh_pred %>%
    mutate(est= case_when(
      trt=="placebo" ~ 1-exp(-hazard*exp(as.vector(pred_data[pred_data$trt == "placebo", "est"]))),
      trt=="cort" ~ 1-exp(-hazard*exp(as.vector(pred_data[pred_data$trt == "cort", "est"])))
    ),
    CI_up= case_when(
      trt=="placebo" ~ 1-exp(-hazard*exp(as.vector(pred_data[pred_data$trt == "placebo", "upperCI"]))),
      trt=="cort" ~ 1-exp(-hazard*exp(as.vector(pred_data[pred_data$trt == "cort", "upperCI"])))
    ),
    CI_low= case_when(
      trt=="placebo" ~ 1-exp(-hazard*exp(as.vector(pred_data[pred_data$trt == "placebo", "lowerCI"]))),
      trt=="cort" ~ 1-exp(-hazard*exp(as.vector(pred_data[pred_data$trt == "cort", "lowerCI"])))
    ))
  
  #remove the point that was added as a trick so that geom_step() displays CIs until last day of treatment
  bh_pred[c(nrow(bh),nrow(bh)*2),'est']<-NA
  
  bh_pred_temp<-bh_pred %>% mutate(date_capt=capdat)
  
  #Create dataframe for figure with all predictions
  if(capdat==selectedcapdates[1]){
    
    bh_pred_fig<-bh_pred_temp
    
  }else{
    
    bh_pred_fig<-bh_pred_fig%>%
      bind_rows(bh_pred_temp)
  }
  
}

#Trick to make CI appear until last point because geom_step() only plots a line until the last point
bh_pred_fig_B<-bh_pred_fig%>%
  filter(date_capt%in%c(7),day_since_trt%in%c(7:19)) #select capture date (cap date 7 = May 6th) and period of interest to display
bh_pred_fig_B[c(nrow(bh_pred_fig_B)/2,nrow(bh_pred_fig_B)),'est']<-NA

# plot
library(ggridges)
B <-ggplot(bh_pred_fig_B , aes(day_since_trt-0.5, est, color= trt)) +
  geom_step()+
  geom_point(aes(x=day_since_trt, est, color= trt))+
  pammtools::geom_stepribbon(aes(ymin = CI_low, ymax = CI_up, fill = trt), alpha = .3)+
  labs(y= "Migratory departure probability", x = "Day after pellet implantation", fill="", color="")+
  scale_x_continuous(n.breaks=12)+
  scale_color_manual(values=c("#E73000","#529CB5"),labels=c('CORT','PLACEBO'))+
  scale_fill_manual(values=c("#E73000","#529CB5"),labels=c('CORT','PLACEBO'),
                    guide = guide_legend(override.aes = list(color = NA)))+
  theme_classic()+
  theme(text = element_text(size = 20),
        legend.position='top')
B


#.....................................
## Effect size for May 6 (date_capt=7)
#.....................................
bh_pred_fig%>%filter(date_capt%in%c(7),day_since_trt%in%c(14))

#.....................................
## Figure 2. Combined A and B
#.....................................
library(gridExtra)
A_B <- grid.arrange(A, B, nrow = 1)
ggsave("Figure_departure.png", A_B, width = 15, height = 7, dpi=300)
ggsave("Figure_2.pdf", A_B, width = 15, height = 7)

#......................
## Appendix. Figure S6.
#......................

# Add real dates labels to dataframe
first_capture_yday <- unique(instance_df[instance_df$date_capt==1, "date_deploiement"])
first_capture_date <- as.Date("2021-01-01")+first_capture_yday
bh_pred_fig$date_capt_long = first_capture_date+as.numeric(as.character(bh_pred_fig$date_capt))-1
bh_pred_fig$current_date_yday = first_capture_yday+as.numeric(as.character(bh_pred_fig$date_capt))+bh_pred_fig$day_since_trt

# relevel capture so that they are nicely displayed
l=rev(c(1,2,4,5,7,9,10,11,13,14))
bh_pred_fig$date_capt <- factor(bh_pred_fig$date_capt, levels = l)

# Set date labels
bh_pred_label <- unique(bh_pred_fig[, c("date_capt", "date_capt_long")])
Sys.setlocale("LC_TIME", "en_US.UTF-8")
bh_pred_label$capt_labels = format(bh_pred_label$date_capt_long, "Treatment on %B %d")
bh_pred_fig$current_date <- as.Date("2021-01-01")+bh_pred_fig$current_date_yday




# plot
C <- ggplot(bh_pred_fig, aes(x=current_date_yday-0.5, est, color= trt)) +
  geom_step()+
  geom_point(aes(x=current_date_yday))+
  pammtools::geom_stepribbon(aes(ymin = CI_up, ymax = CI_low, fill = trt), alpha = .3)+
  geom_text(aes(x = -Inf, y = 0.9, hjust = -0.1, label = capt_labels), size=4.5, color="black", data = bh_pred_label, vjust = 1) +
  
  scale_color_manual(values=c("#e9451a","#579db9"), labels=c('Corticosterone','Placebo'))+
  scale_fill_manual(values = c("#e9451a","#579db9"), labels=c('Corticosterone','Placebo'))+
  
  scale_x_continuous(n.breaks = 11, expand = c(0, 0))+
  #scale_x_continuous(breaks = seq(2,22, 2), expand = c(0, 0))+
  scale_y_continuous(breaks = c(0.25, 0.5, 0.75, 1))+
  geom_hline(yintercept=0)+
  
  coord_cartesian(xlim=c(132,144.5))+
  labs(y= "Departure Probability", x = "Day of year", fill="", color="")+
  theme_classic(base_size = 10)+
  theme(text = element_text(size = 20),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=14),
        axis.title.y=element_text(margin =margin(r=5)),
        legend.position = "top")+
  facet_wrap(date_capt~., ncol=1)+
  theme(strip.text = element_blank())
C

ggsave("Figure_annexe_S6.png", C, width = 7.8, height = 10.5, dpi=400)



# TABLE 2. 
library(coxme)
instance_df <- instance_df %>% mutate(trt=factor(trt, levels = c("placebo", "cort")))## Relevel factors to obtain  CORT estimates

summary(cme)
  
p <- sjPlot::tab_model(cme, dv.labels = "Probability of departure",transform=NULL)
p

#..................................................
### Computation of R2 for departure date model ####
#..................................................

#Compute some approximation of R2 based on likelihood ratio tests
library(CoxR2)

#modify coxr2 function from CoxR2 package to work with cox mixed effects objects
coxmer2 <- function (object){
  cox <- object
  
  #################################### #
  rval <- list()
  #################################### #
  if (!is.null(cox$n[1])) #This was modified from original function: if (!is.null(cox$nevents))
    rval$nevent <- cox$n[1] #This was modified from original function: rval$nevent <- cox$nevent
  
  df <-   cox$df[1] #For the actual df and not those at the penalized location
  ##################################################### #
  logtest <- -2 * (cox$loglik[1] - cox$loglik[2]) #Difference between null model likelihood and integrated partial likelihood at the solution
  ##################################################### #
  rval$logtest <- c(test = logtest, df = df, pvalue = pchisq(logtest,
                                                             df, lower.tail = FALSE))
  
  ##################################################################################### #
  rval$rsq <- c(rsq = 1 - exp(-logtest/cox$n[1]))  ###n <- number of uncensored events
  ###################################################################################### #
  
  rval
}
#run model without random effects
cph <- coxph(Surv(day_since_trt, is_departed) ~ trt*date_capt, data= instance_df, x=T)

#Compute 'R2' statistics
(R2_marginal <- coxr2(cph)$rsq)       #'r2' of model without random effect, i.e. something similar to a 'marginal R2'= 0.92
(R2_conditional <- coxmer2(cme)$rsq)  #'r2' of model based on LRT, something analogous to a 'conditional R2'= 0.95

#Note that the authors of the package have removed these statistics from the summary because of their unreliability.




# —————————————————————————————————————————————————————————————————————————————
# Analysis 3. Comparison of habitat use by cort- and placebo-treated geese ####
# —————————————————————————————————————————————————————————————————————————————

library(sf)
library(raster)
library(mapview)

## Format data
gps_habitat_data<-readRDS("gps_landclass.RDS")

# drop geometry now that we don't need it
exp_std_all<-st_drop_geometry(gps_habitat_data)

#Removing pair number 9 because HJ was injured during the capture
exp_std_all <- exp_std_all %>% filter(paire!=9)

#Filter all point that are in the St Laurent region
exp_std_stlo<-exp_std_all%>%
  filter(in_stlo=='IN')

#This restricts data for each pair to days when both members of the pair are in the St-Lawrence valley
exp_std_restr<-exp_std_all%>%
  filter(jday<=last_stlo,
         jday>=start_stlo)

#Select which dataset to analyse
dataset<-exp_std_restr

#St-Lawrence data
#Compute number of points taken every day
fixperday<- dataset%>%
  group_by(ID, jday)%>%
  dplyr::summarize(ntot=n())%>%ungroup()

#Compute proportion of points that are in each habitat 
prop_time<-dataset%>%
  mutate(TOD=as.numeric(hour(local_timestamp)))%>%
  mutate(periode=case_when(TOD<8~'MORNING',
                           TOD>17~'EVENING',
                           TOD%in%c(8:17)~'DAY'),
         crop_num=ifelse(class_agri=='cropland',1,0))%>%
  left_join(fixperday, by=c('ID','jday'))%>%
  group_by(ID, jday, JT, trt, paire, periode)%>% #Grouping by class allows computing the number of points in each habitat class
  dplyr::summarize(total_points = n(),
                   pct=sum(crop_num)/total_points)

#Keep only days of treatment 2:10 (for which we confirmed an effect of the CORT implant in captivity)
prop_crop<-prop_time%>%
  filter(JT<=10)

prop_crop$trt <- factor(prop_crop$trt, levels=c('placebo', 'cort'))
prop_crop$periode <- factor(prop_crop$periode, levels=c('DAY', 'EVENING', 'MORNING'))


#.............................................................
### Generalized linear mixed model (GLMM) for habitat use ####
#.............................................................

## Model fitting ####
mod_habitat<-MASS::glmmPQL(pct~trt * JT + periode, random = list(~1|paire, ~1|ID), family='quasibinomial', weights=total_points, data=prop_crop[prop_crop$JT>1,])
summary(mod_habitat)

#Get confidence intervals for habitat use parameters
I=intervals(mod_habitat, which="fixed")
d=data.frame(I$fixed)
round(d[,c('est.','lower','upper')],2)

#...................................
### Predictions and effect size ####
#...................................

# build newdata for predictions
newdata <- data.frame(expand.grid(pct=1, 
                                  trt = levels(prop_crop$trt),
                                  periode= levels(prop_crop$periode),
                                  JT = seq(from= min(prop_crop$JT), to= max(prop_crop$JT), by=1)))

newdata <- newdata %>%
  mutate(trt=factor(newdata$trt, levels=c("placebo", "cort"))) %>% 
  arrange(trt, JT)

# model matrix for newdata
X_new <- model.matrix(formula(mod_habitat), data = newdata)

# fixed effect estimates
beta <- fixef(mod_habitat)

# logit predictions
newdata$pred_logit <- X_new %*% beta

# standard errors
vcov_mat <- vcov(mod_habitat) # Extract variance-covariance matrix
se_logit <- sqrt(diag(X_new %*% vcov_mat %*% t(X_new))) 

# 95% confidence intervals on logit scale
newdata$logit_low <- newdata$pred_logit - 1.96 * se_logit
newdata$logit_high <- newdata$pred_logit + 1.96 * se_logit

# probability scale
newdata$pred_prob <- plogis(newdata$pred_logit)
newdata$prob_low <- plogis(newdata$logit_low)
newdata$prob_high <- plogis(newdata$logit_high)

# filter to keep only day for graph
newdata_plot <- newdata[newdata$periode=='DAY' & newdata$JT>=2,] 

#.......................................................
### Fig. Foraging intensity of greater snow goose ####
#.......................................................

prop_crop_sum <- prop_crop %>% 
  filter(periode=="DAY") %>% 
  group_by(JT, trt) %>% 
  summarize(prop_cropland = mean(pct),
            se= sd(pct)/sqrt(n())) %>% 
  mutate(prop_low=prop_cropland-1.96*se,
         prop_high=prop_cropland+1.96*se
  )

dodge_width <- 0.2
newdata_plot$trt <- factor(newdata_plot$trt, levels=c("cort", "placebo")) 
prop_crop_day <- prop_crop[prop_crop$periode == "DAY",]
prop_crop_day$trt <- factor(prop_crop_day$trt, levels=c("cort", "placebo"))
prop_crop_day$group <- ifelse(prop_crop_day$trt=="cort", paste0("B",prop_crop_day$JT), paste0("A",prop_crop_day$JT))

Z <- ggplot(prop_crop_day, aes(x= pct, y= JT-dodge_width, fill = trt)) +
  
  geom_path(data=newdata_plot, aes(y=JT-dodge_width, x=prob_low, col=trt), show.legend=F, lty=2, alpha=1)+
  geom_path(data=newdata_plot, aes(y=JT-dodge_width, x=prob_high, col=trt), show.legend=F, lty=2, alpha=1)+
  
  # display lines
  geom_path(data=newdata_plot, aes(y=JT-dodge_width, x=pred_prob, col=trt), show.legend=F, lty=1, lwd=1.2, alpha=1)+
  
  # display point and errorbar
  geom_errorbar(data=prop_crop_sum, aes(y=JT-dodge_width, x=prop_cropland, xmin= prop_low, xmax=prop_high),
                show.legend = F, alpha=1, width=0, size=1)+
  
  geom_point(data = subset(prop_crop_sum, trt == "placebo"), aes(y = JT - dodge_width, x = prop_cropland, fill = trt),
             show.legend = F, alpha = 1, size = 4.5, stroke = 1.7, pch = 21) +
  
  geom_point(data = subset(prop_crop_sum, trt == "cort"), aes(y = JT - dodge_width, x = prop_cropland, fill = trt),
             show.legend = F, alpha = 1, size = 4.5, stroke = 1.7, pch = 21) +
  scale_fill_manual(values = c("#F3977F","#A9CDDB"), breaks = c("cort", "placebo"), labels= c("CORT", "placebo"), guide = "none")+
  
  # new scale to display violin plot
  new_scale_fill() +
  geom_density_ridges2(data=prop_crop_day, aes(y= JT, x= pct, fill= trt, group=group),
                       scale =0.7, show.legend = T, alpha=0.4, bandwidth = 0.045, col= NA)+
  
  scale_y_continuous(n.breaks=10)+
  coord_flip(xlim= c(0,1))+
  
  xlab("Proportion of time in fields")+
  ylab("Days after pellet implantation")+
  theme_minimal()+ 

  scale_fill_manual(values=c("#E73000", "#529CB5"), breaks = c("cort", "placebo"), labels= c("CORT", "PLACEBO"),
                    guide = guide_legend(override.aes = list(color = NA)))+
  scale_color_manual(values = c("#F3977F","#A9CDDB"), breaks = c("cort", "placebo"), labels= c("CORT", "PLACEBO"))+
  
  theme(panel.grid=element_blank(),
        #panel.background = element_rect(fill='white'), #transparent panel bg
        #plot.background = element_rect(fill='white'),
        plot.title = element_text(size=22, face="bold", hjust = -0.3),
        legend.title=element_blank(),
        legend.text = element_text(size=14, colour="black"),
        axis.text=element_text(size=17, colour="black"),
        axis.title=element_text(margin =margin(0,10,20,), size=18,face="bold"),
        axis.line = element_line(color = 'black'),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 8)))
Z

ggsave("Figure_cropland.png", Z, width = 10, height = 5.8, dpi=300, bg = 'white')


# —————————————————————————————————————————
# FIGURE S2: quantitative data on CORT ####
# —————————————————————————————————————————

# graylag
graylag <- readRDS("graylag_CORT.RDS")
graylag$Corticostérone.ng.mL <- stringr::str_replace(graylag$Corticostérone.ng.mL, ',', '.')
graylag$CORT <- as.numeric(graylag$Corticostérone.ng.mL)
graylag_trt <- graylag[graylag$TRT=='CORT' & graylag$JT>=0 & graylag$JT!=4,]

graylag_baseline <- graylag[graylag$TRT=='BASELINE',]

lower_BL <- as.numeric(quantile(graylag_baseline$CORT, c(0.025,0.975))[1])
upper_BL <- as.numeric(quantile(graylag_baseline$CORT, c(0.025,0.975))[2])
mean_BL <- mean(graylag_baseline$CORT)

# sngo
sngo <- readRDS("sngo_CORT.RDS")
sngo$basal_CORT <- as.numeric(sngo$basal_CORT)
lower_BL_sngo <- as.numeric(quantile(sngo$basal_CORT, c(0.025,0.975), na.rm=T)[1])
upper_BL_sngo <- as.numeric(quantile(sngo$basal_CORT, c(0.025,0.975), na.rm=T)[2])
mean_BL_sngo <- mean(sngo$basal_CORT, na.rm=T)

sngo$stress_induced_CORT <- as.numeric(sngo$stress_induced_CORT)
lower_STRESS_sngo <- as.numeric(quantile(sngo$stress_induced_CORT, c(0.025,0.975), na.rm=T)[1])
upper_STRESS_sngo <- as.numeric(quantile(sngo$stress_induced_CORT, c(0.025,0.975), na.rm=T)[2])
mean_STRESS_sngo <- mean(sngo$stress_induced_CORT, na.rm=T)



graylag_trt_sum <- graylag_trt %>% 
                    group_by(JT) %>% 
                    summarize(avg_cort=mean(CORT, na.rm=T),
                              se_cort=sd(CORT, na.rm=T)/sqrt(n())) %>% 
                    mutate(lwr_se=avg_cort-se_cort,
                           upr_se=avg_cort+se_cort,
                           lwr_CI=avg_cort-2*se_cort,
                           upr_CI=avg_cort+2*se_cort)

library(ggplot2)
line_width=1.2
col_sngo_BL <- "#532572"
col_graylag_BL <- "#7796A6"
col_sngo_stress <- "#FFA333"

# Create a data frame for the rectangles with labels
rect_data <- data.frame(
  category = c("Graylag goose \nbaseline", "Snow goose \nbaseline", "Snow goose \nstress"),
  ymin = c(lower_BL, lower_BL_sngo, lower_STRESS_sngo),
  ymax = c(upper_BL, upper_BL_sngo, upper_STRESS_sngo),
  fill_color = c(col_graylag_BL, col_sngo_BL, col_sngo_stress)
)

rect_data$category <- factor(rect_data$category, levels=c("Snow goose \nstress", "Graylag goose \nbaseline", "Snow goose \nbaseline"))

cort_lvl_plot <- ggplot()+

        geom_rect(data = rect_data, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = category), 
                  alpha = 0.1) +
  
        labs(x="Days of treatment",
             y="Corticosterone (ng/ml)")+
  
  
         geom_hline(yintercept=mean_BL, lty=1, col=col_graylag_BL, lwd=line_width)+
         geom_hline(yintercept=lower_BL, lty=3, col=col_graylag_BL, lwd=line_width/1.1)+
         geom_hline(yintercept=upper_BL, lty=3, col=col_graylag_BL, lwd=line_width/1.1)+
        
         geom_hline(yintercept=mean_BL_sngo, lty=1, col=col_sngo_BL, lwd=line_width)+
         geom_hline(yintercept=lower_BL_sngo, lty=3, col=col_sngo_BL,lwd=line_width/1.1)+
         geom_hline(yintercept=upper_BL_sngo, lty=3, col=col_sngo_BL, lwd=line_width/1.1)+
        
         geom_hline(yintercept=mean_STRESS_sngo, lty=1, col=col_sngo_stress, lwd=line_width)+
         geom_hline(yintercept=lower_STRESS_sngo, lty=3, col=col_sngo_stress, lwd=line_width/1.1)+
         geom_hline(yintercept=upper_STRESS_sngo, lty=3, col=col_sngo_stress, lwd=line_width/1.1)+
  
  
         geom_errorbar(data=graylag_trt_sum, aes(x=JT, y=avg_cort, ymin=lwr_CI, ymax=upr_CI), width=0)+
         geom_point(data=graylag_trt_sum, aes(x=JT, y=avg_cort,col='CORT_graylag'), size=3.5)+
         geom_line(data=graylag_trt_sum, aes(x=JT, y=avg_cort, col='CORT_graylag'), lwd=0.8)+
  
          scale_fill_manual(name = "Corticosterone Levels", 
                            values = setNames(rect_data$fill_color, rect_data$category)) +
  
          scale_color_manual(values = c("CORT_graylag" = "black"),  
                             labels = c("CORT_graylag" = "CORT-implanted\nGreylag goose"))+
  
          guides(fill = guide_legend(override.aes = list(alpha = 1), byrow=TRUE, order=1),
                 color=guide_legend(order=2))+

    
         scale_x_continuous(n.breaks=10)+
         theme(panel.grid=element_blank(),
               panel.background = element_rect(fill='white'), #transparent panel bg
               plot.background = element_rect(fill='white'),
               plot.title = element_text(size=22, face="bold", hjust = -0.3),
               legend.title=element_blank(),
               legend.text = element_text(size=12, colour="black"),
               legend.background = element_rect(fill="white"),
               legend.position="right",
               legend.key.spacing.y = unit(1,"line"),
               #legend.key.size=unit(2,"line"),
               axis.text=element_text(size=17, colour="black"),
               axis.title=element_text(margin =margin(0,10,20,), size=18,face="bold"),
               axis.line = element_line(color = 'black'),
               axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 8)))

cort_lvl_plot

ggsave("Figure_sup_map_CORT_quantitative_CI.png", cort_lvl_plot, width=10, height=8, dpi=300)








