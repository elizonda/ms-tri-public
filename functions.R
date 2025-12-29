#Functions for Triangle Paper (2025)

# Detection ----

dprime_env_t.test = function(dprime_df, condition = "env") {
  library(tidyverse)
  library(effectsize)
  
  data <- dprime_df
  
  if (condition == "env"){
    
    {char_length1 <- nchar("Remember: Here it shows directionality between the Flat and Triangle tones for dprime scores.")
      cat(red(paste0("\n", strrep("=", char_length1+3), "\n")))
      cat(paste0(white$bold('T-test between Envlopes Comparison: ', green$bold$underline('Triangle against Flat', "\n"))))
      cat(blue$italic("Remember: Here it shows directionality between the Flat and Triangle tones for dprime scores. \n"))
      cat(blue$italic("If significant, it means that the triangle is better than the flat tone for detectability!"))
      cat(red(paste0("\n", strrep("=", char_length1+3), "\n")))
    } #Character print options to look pretty :)
    
    data %>% 
      dplyr::filter(env =="PercTri") -> triangle
    
    data %>% 
      dplyr::filter(env =="Flat") -> flat
    
    cat(green("T-Test"))
    print(t.test(triangle$dprime, flat$dprime, alternative = "greater", paired = T))
    
    cat(green("Paired Cohen's D \n"))
    percVflat <- triangle$dprime - flat$dprime
    print(abs(mean(percVflat) / sd(percVflat)))
    
    cat(green("Degrees of Freedom (Paired) \n"))
    print((length(temp$diff)/2) - 1)
  }
  else {
    data %>% 
      dplyr::filter(env =="PercTri")  %>% 
      dplyr::group_by(participant) %>% 
      dplyr::summarize(
        diffs
      ) -> triangle
    
    
  }
}

dprime_MixedEffects_ANOVA = function(dprime_df1, dprime_df2, dprime_df3) {
  library(afex)
  library(emmeans)
  library(stringr)
  library(effectsize)
  
  columns_to_keep <- c("participant", "snr", "env", "dprime")
  
  temp_df1 <- dprime_df1 %>% 
    dplyr::select(all_of(columns_to_keep)) %>% 
    mutate(experiment = "1a",
           participant = str_c(participant, "_1a"))
  
  temp_df2 <- dprime_df2 %>% 
    dplyr::select(all_of(columns_to_keep)) %>% 
    mutate(experiment = "1b",
           participant = str_c(participant, "_1b"))
  
  temp_df3 <- dprime_df3 %>% 
    dplyr::select(all_of(columns_to_keep)) %>% 
    mutate(experiment = "1c",
           participant = str_c(participant, "_1c"))
  
  
  # Combine
  combined_df <- rbind(temp_df1, temp_df2, temp_df3)
  
  mixed_model <- aov_car(
    dprime ~ env * snr * experiment + Error(participant/(env * snr)), 
    data = combined_df
  )
  
  print(anova(mixed_model))
  
  # Test interactions (e.g., does condition effect vary by experiment?)
  cat(green("Pairwise ~ env | Experiment \n"))
  print(emmeans(mixed_model, pairwise ~ env | experiment))
  
  #Effect sizes
  print(effectsize::eta_squared(mixed_model, partial = T))
  
  # Plot marginal means
  print(emmip(mixed_model, env ~ experiment, CIs = TRUE) +
          theme_classic())
  
  #Tukey Pairs
  emm_interaction <- emmeans(mixed_model, ~ env * experiment)
  
  tukey_exp_within_env <- pairs(
    emm_interaction,
    simple = "experiment",  # Compare experiments within each env level
    adjust = "tukey"
  )
  
  print(summary(tukey_exp_within_env, infer = TRUE))
}

dprime_effectsize_META_test_between = function(dprime_df1, dprime_df2){
  library(metafor)
  library(crayon)
  
  n1 <- length(unique(dprime_df1$participant))
  
  n2 <- length(unique(dprime_df2$participant))
  
  aov_df1 <- aov_car(dprime~snr*env + Error(participant/(snr*env)), 
                         data=dprime_df1)
  
  aov_df2 <- aov_car(dprime~snr*env + Error(participant/(snr*env)), 
                     data=dprime_df2)
  
  eta_sq1 <- effectsize::eta_squared(aov_df1, partial = T, alternative = "two.sided")$Eta2_partial[2]
  eta_sq2 <- effectsize::eta_squared(aov_df2, partial = T, alternative = "two.sided")$Eta2_partial[2]
  
  # Example: η²ₚ from two independent ANOVAs
  study1 <- list(eta_sq = eta_sq1, n = n1)  # ANOVA 1
  study2 <- list(eta_sq = eta_sq2, n = n2)  # ANOVA 2

  # Convert η²ₚ to Cohen's f and SE
  data <- data.frame(
    study = c("study 1", "study 2"),
    f = sqrt((c(study1$eta_sq, study2$eta_sq)) / (1 - c(study1$eta_sq, study2$eta_sq))),
    se = 1 / (c(study1$n, study2$n)))
  
  
  # Meta-analysis (random-effects model)
  meta_model <- rma(yi = f, sei = se, data = data, method = "REML")
  
  # Test for heterogeneity (Q-test)
  print(summary(meta_model))
  
  cat(green(underline("Remember:\n")))
  cat(silver("Q statistic determines heterogeneity -> if not singificant then", red(bold("η²ₚ is consistent across experiments.\n"))))
  cat(silver("τ² statistic determines true variance -> if greater than 0 then", red(bold("heterogeneity exists between η²ₚ experiments. \n"))))
  cat(silver("I² what percentage of total variability is due to true differences: \n", 
             red(bold("0%: All variability is from sampling error \n", "25%/50%/75%: Low/medium/high heterogeneity."))))
}


# Speech Comprehension ----

## Stats ----

speech_MeanSD = function(data) {
  library(tidyverse)
  library(afex)
  library(crayon)
  
  speech_temp_AggData <- data  %>% 
    dplyr::group_by(participant, env, snr) %>% 
    dplyr::summarize(
      CRMhit = mean(CRMhit)
    ) %>% 
    dplyr::ungroup() 
  
  speech_temp_MeanSD_SnrEnv <- speech_temp_AggData %>%
    dplyr::group_by(env, snr) %>% 
    dplyr::summarize(mean_hit = mean(CRMhit),
                     sd_hit = sd(CRMhit))
  
  speech_temp_MeanSD_env <- speech_temp_AggData %>%
    dplyr::group_by(env) %>% 
    dplyr::summarize(mean_hit = mean(CRMhit),
                     sd_hit = sd(CRMhit))
  
  speech_temp_MeanSD_snr <- speech_temp_AggData %>%
    dplyr::group_by(snr) %>% 
    dplyr::summarize(mean_hit = mean(CRMhit),
                     sd_hit = sd(CRMhit))
  
  
  {char_length1 <- nchar("Means & SDs of Signal-to-Noise Ratio vs Envelope (tone type)")
    cat(red(paste0("\n", strrep("=", char_length1+3), "\n")))
    cat(paste0(white$bold('Means & SDs of ', green$bold$underline('Signal-to-Noise Ratio vs Envelope (tone type)'))))
    cat(red(paste0("\n", strrep("=", char_length1+3), "\n")))
  } #Character print options to look pretty :)
  
  print(speech_temp_MeanSD_SnrEnv)
  
  {char_length2 <- nchar("Means & SDs of Envelope (tone type)")
    cat(red(paste0("\n", strrep("=", char_length2+3), "\n")))
    cat(paste0(white$bold('Means & SDs of ', green$bold$underline('Envelope (tone type)'))))
    cat(red(paste0("\n", strrep("=", char_length2+3), "\n")))
  } #Character print options to look pretty :)
  
  print(speech_temp_MeanSD_env)
  
  {char_length3 <- nchar("Means & SDs of Signal-to-Noise Ratio")
    cat(red(paste0("\n", strrep("=", char_length3+3), "\n")))
    cat(paste0(white$bold('Means & SDs of ', green$bold$underline('Signal-to-Noise Ratio'))))
    cat(red(paste0("\n", strrep("=", char_length3+3), "\n")))
  } #Character print options to look pretty :)
  
  print(speech_temp_MeanSD_snr)
  
}

speech_rmANOVA = function(data, effect_size = T, env_exclusive = T) {
  library(tidyverse)
  library(afex)
  library(crayon)
  
  speech_temp_AggData <- data  %>% 
    dplyr::group_by(participant, env, snr) %>% 
    dplyr::summarize(
      CRMhit = mean(CRMhit)
    ) %>% 
    dplyr::ungroup() 
  
  speech_temp_aov <- aov_car(CRMhit~snr*env + Error(participant/(snr*env)),
          data=(speech_temp_AggData 
                %>% filter(env != "silence") #use AOV with Silence, use aov_car without Silence
                ))
  
  {char_length1 <- nchar("Repeated Measures ANOVA:  Signal-to-Noise ratio vs. Envelope (tone type)")
    cat(red(paste0("\n", strrep("=", char_length1+3), "\n")))
    cat(paste0(white$bold('Repeated Measures ANOVA: ', green$bold$underline('Signal-to-Noise ratio vs. Envelope (tone type)'))))
    cat(red(paste0("\n", strrep("=", char_length1+3), "\n")))
  } #Character print options to look pretty :)
  
  print(summary(speech_temp_aov))
  
  if(effect_size == T){
    {char_length2 <- nchar("Effect Size:  Signal-to-Noise ratio vs. Envelope (tone type)")
      cat(red(paste0("\n", strrep("=", char_length2+3), "\n")))
      cat(paste0(white$bold('Effect Size: ', green$bold$underline('Signal-to-Noise ratio vs. Envelope (tone type)'))))
      cat(red(paste0("\n", strrep("=", char_length2+3), "\n")))
    } #Character print options to look pretty :)
    print(effectsize::eta_squared(speech_temp_aov, alternative = "two.sided"))
  }
  
  if(env_exclusive == T){
    CRMperformance_tri_only <- speech_temp_AggData %>% 
      dplyr::filter(env == "PercTri")
    
    aov_CRMperformance_tri_only <- aov_car(CRMhit~snr + Error(participant/snr), data = CRMperformance_tri_only)
    
    {char_length3 <- nchar("Repeated Measures ANOVA: Triangle Alone")
      cat(red(paste0("\n", strrep("=", char_length3+3), "\n")))
      cat(paste0(white$bold('Repeated Measures ANOVA: ', green$bold$underline('Triangle Alone'))))
      cat(red(paste0("\n", strrep("=", char_length3+3), "\n")))
    } #Character print options to look pretty :)
    
    print(summary(aov_CRMperformance_tri_only))
    
    if(effect_size == T){
      
      {char_length4 <- nchar("Effect Size: Triangle Alone")
      cat(red(paste0("\n", strrep("=", char_length4+3), "\n")))
      cat(paste0(white$bold('Effect Size: ', green$bold$underline('Triangle Alone'))))
      cat(red(paste0("\n", strrep("=", char_length4+3), "\n")))
      } #Character print options to look pretty :)
      
      print(effectsize::eta_squared(aov_CRMperformance_tri_only, alternative = "two.sided"))
    }
    
    CRMperformance_flat_only <- speech_temp_AggData %>% 
      dplyr::filter(env == "Flat")
    aov_CRMperformance_flat_only <- aov_car(CRMhit~snr + Error(participant/snr), data = CRMperformance_flat_only)
    
    {char_length5 <- nchar("Repeated Measures ANOVA: Flat Alone")
      cat(red(paste0("\n", strrep("=", char_length5+3), "\n")))
      cat(paste0(white$bold('Repeated Measures ANOVA: ', green$bold$underline('Flat Alone'))))
      cat(red(paste0("\n", strrep("=", char_length5+3), "\n")))
    } #Character print options to look pretty :)
    
    print(summary(aov_CRMperformance_flat_only))
    
    if(effect_size == T){
      {char_length6 <- nchar("Effect Size: Flat Alone")
      cat(red(paste0("\n", strrep("=", char_length6+3), "\n")))
      cat(paste0(white$bold('Effect Size: ', green$bold$underline('Flat Alone'))))
      cat(red(paste0("\n", strrep("=", char_length6+3), "\n")))
      } #Character print options to look pretty :)
      
      print(effectsize::eta_squared(aov_CRMperformance_flat_only, alternative = "two.sided"))
    }
  }
}

speech_baseline_t.test = function(data, condition = "env", plot = T) {
  library(tidyverse)
  library(effectsize)
  
  data %>% 
    dplyr::filter(env == "silence") %>% 
    dplyr::group_by(participant) %>% 
    dplyr::summarize(
      silence_hit_rate = (sum(CRMhit)/108)*100,
    ) -> silence_speechRate
  
  tone_baseline <- data %>% 
    dplyr::filter(env != "silence") %>% 
    dplyr::group_by(participant, snr, env, SnrEnv) %>% 
    dplyr::summarize(
      tone_hit_rate = (sum(CRMhit)/9)*100,
    ) %>% 
    left_join(silence_speechRate, by = "participant") %>% 
    dplyr::mutate(diff = (tone_hit_rate - silence_hit_rate)) %>% 
    ungroup()
  
  
  if (condition == "env"){
      
      {char_length1 <- nchar("Remember: if 95% CI contains 0, then the triangle tone is no different than baseline (Silence)")
      cat(red(paste0("\n", strrep("=", char_length1+3), "\n")))
      cat(paste0(white$bold('T-test Baseline Comparison: ', green$bold$underline('Triangle vs. Baseline (Silence)', "\n"))))
      cat(blue$italic("Remember: if 95% CI contains 0, then the triangle tone is no different than baseline (Silence) \nMeaning the triangle tone did not harm speech comprehension"))
      cat(red(paste0("\n", strrep("=", char_length1+3), "\n")))
      } #Character print options to look pretty :)
      
      tone_baseline %>% 
        dplyr::filter(env =="PercTri") -> temp
      cat(green("T-Test"))
      print(t.test(temp$diff, alternative = "less"))

      cat(green("One Sample Cohen's D \n"))
      print(((mean(temp$diff - 0))/sd(temp$diff))) #Single Cohen's D 
      
      cat(green("Degrees of Freedom (Paired) \n"))
      print((length(temp$diff)/2) - 1)
    
    
    
      {char_length2 <- nchar("Remember: if 95% CI contains 0, then the flat tone is no different than baseline (Silence)")
      cat(red(paste0("\n", strrep("=", char_length2+3), "\n")))
      cat(paste0(white$bold('T-test Baseline Comparison: ', green$bold$underline('Flat vs. Baseline (Silence)', "\n"))))
      cat(blue$italic("Remember: if 95% CI contains 0, then the flat tone is no different than baseline (Silence) \nMeaning the flat tone did not harm speech comprehension"))
      cat(red(paste0("\n", strrep("=", char_length2+3), "\n")))
      } #Character print options to look pretty :)
      
      tone_baseline %>% 
        dplyr::filter(env =="Flat") -> temp
      cat(green("T-Test"))
      print(t.test(temp$diff, alternative = "less"))

      cat(green("One Sample Cohen's D \n"))
      print(((mean(temp$diff - 0))/sd(temp$diff))) #Single Cohen's D 
      
      cat(green("Degrees of Freedom (Paired) \n"))
      print((length(temp$diff)/2) - 1)
      
      {char_length3 <- nchar("Remember: Here it shows directionality between the Flat and Triangle tone")
        cat(red(paste0("\n", strrep("=", char_length3+3), "\n")))
        cat(paste0(white$bold('T-test between Envlopes Comparison: ', green$bold$underline('Flat vs. Triangle', "\n"))))
        cat(blue$italic("Remember: Here it shows directionality between the Flat and Triangle tone"))
        cat(red(paste0("\n", strrep("=", char_length3+3), "\n")))
      } #Character print options to look pretty :)
      
      tone_baseline %>% 
        dplyr::filter(env =="PercTri") -> temp1
      
      tone_baseline %>% 
        dplyr::filter(env =="Flat") -> temp2
      
      cat(green("T-Test"))
      print(t.test(temp1$diff, temp2$diff, alternative = "less", paired = T))
      
      cat(green("Paired Cohen's D \n"))
      percVflat <- temp1$diff - temp2$diff
      print(abs(mean(percVflat) / sd(percVflat)))
      
      cat(green("Degrees of Freedom (Paired) \n"))
      print((length(temp$diff)/2) - 1)
  }
  
  if (condition == "snr"){
    for (i in 1:length(unique(subset(tone_baseline, env != "silence")$snr))){
      snr_values <- as.character(unique(subset(tone_baseline, env != "silence")$snr))
      snr_levels <- c(6,5,4,3,2,1)
      
      temp <- subset(tone_baseline, snr == snr_values[i])
      
      cat(red$underline(paste0("SNR Value: ", snr_values[i], "\n")))
      cat(red$underline(paste0("SNR Level: ", snr_levels[i], "\n")))
      cat(green("T-Test"))
      print(t.test(temp$diff, alternative = "less"))
      
      cat(green("Single Sample Cohen's D \n"))
      print(((mean(temp$diff - 0))/sd(temp$diff))) #Single Cohen's D 
      
      cat(green("Degrees of Freedom (Paired) \n"))
      print((length(temp$diff)/2) - 1)
    }
  }
  
  if (condition[1] == "SnrEnv"){
    
    for (i in 1:length(unique(subset(tone_baseline, env != "silence")$SnrEnv))){
      SnrEnv_values <- as.character(unique(subset(tone_baseline, env != "silence")$SnrEnv))
      
      temp <- subset(tone_baseline, SnrEnv == SnrEnv_values[i])
      
      cat(red$underline(paste0("SNR x ENV Value: ", SnrEnv_values[i], "\n")))
      cat(green("T-Test"))
      print(t.test(temp$diff, alternative = "less"))
      
      cat(green("Single Sample Cohen's D \n"))
      print(((mean(temp$diff - 0))/sd(temp$diff))) #Single Cohen's D 
      
      cat(green("Degrees of Freedom (Paired) \n"))
      print((length(temp$diff)/2) - 1)
    }
  }
}

speech_effectsize_META_test_between = function(data1, data2){
  library(metafor)
  library(crayon)
  
  speech_temp_AggData1 <- data1  %>% 
    dplyr::group_by(participant, env, snr) %>% 
    dplyr::summarize(
      CRMhit = mean(CRMhit)
    ) %>% 
    dplyr::ungroup() 
  
  speech_temp_AggData2 <- data2  %>% 
    dplyr::group_by(participant, env, snr) %>% 
    dplyr::summarize(
      CRMhit = mean(CRMhit)
    ) %>% 
    dplyr::ungroup() 
  
  n1 <- length(unique(speech_temp_AggData1$participant))
  
  n2 <- length(unique(speech_temp_AggData2$participant))
  
  aov_df1 <- aov_car(CRMhit~snr*env + Error(participant/(snr*env)), 
                     data=speech_temp_AggData1 %>% 
                       filter(env != "silence"))
  
  aov_df2 <- aov_car(CRMhit~snr*env + Error(participant/(snr*env)), 
                     data=speech_temp_AggData2 %>% 
                       filter(env != "silence"))
  
  eta_sq1 <- effectsize::eta_squared(aov_df1, partial = T, alternative = "two.sided")$Eta2_partial[2]
  eta_sq2 <- effectsize::eta_squared(aov_df2, partial = T, alternative = "two.sided")$Eta2_partial[2]
  
  # Example: η²ₚ from two independent ANOVAs
  study1 <- list(eta_sq = eta_sq1, n = n1)  # ANOVA 1
  study2 <- list(eta_sq = eta_sq2, n = n2)  # ANOVA 2
  
  # Convert η²ₚ to Cohen's f and SE
  data <- data.frame(
    study = c("study 1", "study 2"),
    f = sqrt((c(study1$eta_sq, study2$eta_sq)) / (1 - c(study1$eta_sq, study2$eta_sq))),
    se = 1 / (c(study1$n, study2$n)))
  
  
  # Meta-analysis (random-effects model)
  meta_model <- rma(yi = f, sei = se, data = data, method = "REML")
  
  # Test for heterogeneity (Q-test)
  print(summary(meta_model))
  
  cat(green(underline("Remember:\n")))
  cat(silver("Q statistic determines heterogeneity -> if not singificant then", red(bold("η²ₚ is consistent across experiments.\n"))))
  cat(silver("τ² statistic determines true variance -> if greater than 0 then", red(bold("heterogeneity exists between η²ₚ experiments. \n"))))
  cat(silver("I² what percentage of total variability is due to true differences: \n", 
             red(bold("0%: All variability is from sampling error \n", "25%/50%/75%: Low/medium/high heterogeneity."))))
}

## Plots ----

