# install.packages("easypackages")
require("easypackages")
libraries("irtoys","mcIRT","PP",
          "CTT","TAM","psych","sirt")
# Function ----
Optimal.N = function(n = 500, nitem = 20, 
                     correl = 0.60, nsim = 50,
                     seed = runif(1, 1, 10), iter = 100,
                     accept.rmsd = 0.05,
                     prop = 0.30, SIM.Number = 1,
                     ALPHA = 0.80, accept.info = 0.20, File.Number = 1){
  for(zz in 1:nsim){
    # tryCatch() is for catching errors
    # If an operation or simulation cycle does not run properly,
    # this line stops R from ending the entire simulation
      tryCatch({ 
    # Generate Person Parameters ----
        cor   =  acos(correl)                         # corresponding angle
        x1    =  rnorm(n, 0, 1.00)                    # fixed given data
        x2    =  rnorm(n, 0, 1.00)                    # new random data
        X     =  cbind(x1, x2)                        # matrix
        Xctr  =  scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)
        
        Id   =  diag(n)                               # identity matrix
        Q    =  qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
        P    =  tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
        x2o  =  (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
        Xc2  =  cbind(Xctr[ , 1], x2o)                # bind to matrix
        Y    =  Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1
      
        x =  Y[ , 2] + (1 / tan(cor)) * Y[ , 1]       # final new vector
        x = x*10                
        (cor(x1, x))                                # Validate correlation is what you want
      # Classifying True scores as above or below the 75th percentile
        class = x > summary(x)[5]
      # Saving out the rank order values (saved as integers with 1:N
      # 1 being the highest ability; N being the lowest)
        r.x = rank(x)
    # Generate Item Parameters and Sums Scores Step -------------------------------------------------------------
      # Generate item discrimination that is uniformly drawn from two separate distributions
      # One is good and one is not good
        bad = runif(round(nitem*prop), 0.10, 0.50)
        good = runif((nitem - (round(nitem*prop))), 1.5, 2.5)
      # Bind the 'good' and 'bad' item discriminations
        discrimination = c(good, bad)
      # Small SD was chosen to simulate random item discriminations but with all items still being 'good' in their functioning
        difficulty = rnorm(nitem, 0, 0.75) # Create item difficulties about 0 with SD of 0.75
      # Bind Item Parameters into a matrix that is compatiable with R-package
        Parameters.2PL = matrix(cbind(discrimination, difficulty, 0), nitem, 3)
      # Simulated scores using item parameters and ability
        scores.2PL = sim.raschtype(x, b = difficulty, fixed.a = discrimination)
      # Generate summed scores
        sums = rowSums(scores.2PL)
      # Calculate estimate-to-true correlation for Sums condition
        s.x = cor(sums, x)
      # Calculate estimate-to-outcome correlation for Sums condition
        s.y = cor(sums, x1)
      # Generate rank order of the summed scores
        r.s = rank(sums)
      # Calculate those above/below the 75th percentile rank in the Summed Scores condition
        s.class = sums > summary(sums)[5]
      # Run a comparison against the 'True' ranks and those estimated by the Summed Scores condition
      # See the FAQ's on IAmJonathanPark.com for information on how False Positive/False Negatives were calculated
        m = matrix(NA, n, 1)
          for(i in 1:n){
            m[i,] = class[i] - s.class[i]
          }
      # Turning the output into proportions as False Positives
        f.pos.sum = 100 * sum(m == -1)/n
      # Turning the output into proportions as False Negatives
        f.neg.sum = 100 * sum(m == 1)/n
      # Save the original dataset for use in estimating scores with other testing frameworks
        OG2pl = scores.2PL
    # Alpha Derivation ----
        # Run an item analysis on the dataset
          alpha = itemAnalysis(scores.2PL)
        # Establish a prior alpha-level so the next function can stop if Alpha begins to decrease
          prior.a = alpha
        # The following while-loop will iteratively assess Cronbach's 'alpha if removed'
          # remove the item that would improve alpha the most if removed until a >= 0.80
          
          # If the loop gets to a point where all items would reduce alpha if removed, the
          # loop will close.
      while(alpha$alpha < ALPHA){
        scores.2PL = scores.2PL[,-which.max(itemAnalysis(scores.2PL)$itemReport$alphaIfDeleted)]
        (alpha = itemAnalysis(scores.2PL))
        
        if(alpha$alpha > prior.a$alpha){
          alpha = itemAnalysis(scores.2PL)
          prior.a = alpha
        }else if(alpha$alpha < prior.a$alpha){
          alpha$alpha = 1
          message("ALPHA ANALYSES REACHED PEAK. STOPPING PROCESS.")
        }
      }
      alpha = itemAnalysis(scores.2PL)
      # Calculate summed scores using alpha adjusted scale
        sum.scores.alpha = rowSums(scores.2PL)
      # Calculate estimate-to-true correlation for Alpha condition
        a.x = cor(sum.scores.alpha, x)
      # Calculate estimate-to-outcome correlation for Alpha condition
        a.y = cor(sum.scores.alpha, x1)
      # Calculate the 75th percentile of the alpha scores
        s.class = sum.scores.alpha > summary(sum.scores.alpha)[5]
      # Calculate the rank order of the simulated participants when scored with the alpha test
        r.a = rank(sum.scores.alpha)
        m = matrix(NA, n, 1)
          for(i in 1:n){
            m[i,] = class[i] - s.class[i]
          }
      # Calculating the false positive rate as a proportion
        f.pos.alpha = 100 * sum(m == -1)/n
      # Calculating the false negative rate as a proportion
        f.neg.alpha = 100 * sum(m == 1)/n
    # Factor Analysis ----
      # Call in original dataset to ensure you are analyzing original dataset  
        scores.2PL = OG2pl
      # Estimate initial set of factor loadings assuming a unidimensional factor
        # structure
        f.scores = factanal(scores.2PL, 1, scores = "Bartlett")
      # The following while-loop will iteratively remove items until all items
        # possess a factor loading above 0.30; obviously this can be adjusted for
        # more or less conservative estimates/tests but this is a generally accepted
        # heuristic.
          while(min(f.scores$loadings) < 0.30){
            scores.2PL = scores.2PL[,-which.min(f.scores$loadings)]
            f.scores = factanal(scores.2PL, 1, scores = "Bartlett")
            # What to do if CFA deletes all items and prevent the code from crashing
              if(ncol(scores.2PL) < 10){
                f.scores$loadings = 0.5
                message("CFA CANNOT SAVE SCALE USING ADJUSTMENT")
        }
      }
      # Be sure that the final estimate carries over from the while-loop
        # This step is probably redundant but I'd rather be safe than sorry
          f.scores = factanal(scores.2PL, 1, scores = "Bartlett")
      # Turn the initial value into a vector of factor-scores  
        f.scores = f.scores$scores
      # Calculate estimate-to-true correlation for CFA condition
        f.x = cor(f.scores, x)
      # Calculate estimate-to-outcome correlation for CFA condition
        f.y = cor(f.scores, x1)
      # Calculate the rank order of participants when scored using factor analysis
        r.f = rank(f.scores)
      # Generate the 75th percentile rank of participants
        s.class = f.scores > summary(f.scores)[5]
        m = matrix(NA, n, 1)
        for(i in 1:n){
          m[i,] = class[i] - s.class[i]
        }
      # Calculate false positive rate  
        f.pos.cfa = 100 * sum(m == -1)/n
      # Calculate false negative rate
        f.neg.cfa = 100 * sum(m == 1)/n
    # IRT Analysis ---------------------------------------------------------------------
      scores.2PL = OG2pl
        options(warn=-1)
      # Conduct an exploratory factor analysis to receive/calculate the
      # the proportion of variance accounted for by the 'general' factor
      # See, Preston, 2018 (Presentation at WPA) for more information
        rfac = suppressMessages(tetrachoric(scores.2PL,na.rm=TRUE)$rho)
        g3 = suppressMessages(schmid(rfac, nfactors = 3, fm = "minres",digits=3,rotate="oblimin"))
        sumout = colSums(g3$sl[,1:(ncol(g3$sl)-3)])^2
        ECV2 = sumout[1] / sum(sumout)
      
      # Eliminate Items with low loadings (e.g., < 0.30) on the general factor
        while(g3$sl[,1][which.min(abs(g3$sl[,1]))] < 0.30){
          scores.2PL = scores.2PL[,-which.min(abs(g3$sl[,1]))]
          rfac = suppressMessages(tetrachoric(scores.2PL,na.rm=TRUE)$rho)
          g3 = suppressMessages(schmid(rfac, nfactors = 3, fm = "minres",digits=3,rotate="oblimin"))
          (sumout = colSums(g3$sl[,1:(ncol(g3$sl)-3)])^2)
          (ECV2 = sumout[1] / sum(sumout))
        }
      # Estimate the 2PL on the original scores
        estimate.2PL = tam.mml.2pl(scores.2PL, verbose = FALSE)
      # Obtain fit estimates
        fit.2PL = IRT.itemfit(estimate.2PL)
      # Obtain the maximum-achieved information for each item
      # This could have been done in far fewer steps...
        item.information = IRT.informationCurves(estimate.2PL)
        item.information = item.information$info_curves_item
        item.information = rowMaxs.sirt(item.information)$maxval
      # The following while loops comprise two steps.
        # 1. Remove items with the lowest test information.
          # This step is repeated iteratively until:
            # All items have information above 0.20
            # The Log-Likelihood test becomes statistically significant
            # The loop deletes all items (Unlikely)
              a = 1           
              Dif.Test = 0
              while(a == 1){
              # Look at lowest item info for removal
                while(min(item.information) < accept.info){
                  # Generate a model without the worst item (defined as lowest max info.)
                    scores.2PL = scores.2PL[,-which.min(item.information)]
                    mod.fit = estimate.2PL$ic
                  # Save characteristics of the original model (i.e., Log-Likelihood and DF)
                    old.estimate = estimate.2PL; LL1 = mod.fit[3]; NPars1 = mod.fit[11]
                  # Calculate the model characteristics of this new model
                    estimate.2PL = tam.mml.2pl(scores.2PL, verbose = FALSE)
                    fit.2PL = IRT.itemfit(estimate.2PL)
                  # Save item information for the new model  
                    item.information = IRT.informationCurves(estimate.2PL)
                    item.information = item.information$info_curves_item
                    item.information = rowMaxs.sirt(item.information)$maxval
                  # Compare the two nested models
                    mod.fit = estimate.2PL$ic; LL2 = mod.fit[3]; NPars2 = mod.fit[11]
                    LR.test = as.numeric(2*(LL1 - LL2))
                    DF = as.numeric(NPars1 - NPars2)
                    Dif.Test = pchisq(LR.test, DF, lower.tail = TRUE)
                }
              # Remove items with high RMSD (> 0.05)
              # This generally ensures that item expected characteristics closely
                # resemble the observed characteristics. Departures from this may
                # be due to DIF.
                  while(max(fit.2PL$rmsea) > accept.rmsd){
                    # Generate a model without the worst item (defined as highest RMSD)
                      scores.2PL = scores.2PL[,-which.max(fit.2PL$rmsea)]
                      mod.fit = estimate.2PL$ic
                    # Save characteristics of the original model (i.e., Log-Likelihood and DF)
                      old.estimate = estimate.2PL; LL1 = mod.fit[3]; NPars1 = mod.fit[11]
                    # Calculate the model characteristics of this new model
                      estimate.2PL = tam.mml.2pl(scores.2PL, verbose = FALSE)
                      fit.2PL = IRT.itemfit(estimate.2PL)
                    # Compare the two nested models
                      mod.fit = estimate.2PL$ic; LL2 = mod.fit[3]; NPars2 = mod.fit[11]
                      LR.test = as.numeric(2*(LL1 - LL2))
                      DF = as.numeric(NPars1 - NPars2)
                      Dif.Test = pchisq(LR.test, DF, lower.tail = TRUE)
                  }
              # This line forces the end of the loop
                if(min(item.information) > accept.info && max(fit.2PL$rmsea) < accept.rmsd){a = 2}
      }
              options(warn=0)
              estimate.2PL = tam.mml.2pl(scores.2PL, verbose = FALSE)
      # Calculate the final test information
        item.information = IRT.informationCurves(estimate.2PL)
      # Create a plot of the final test information
        plot(item.information$theta, (item.information$test_info_curve)/ncol(scores.2PL),
           type = "l", ylim = c(0,1),
           xlab = "Theta", ylab = "Information", main = "IRT Scores Information")
      # Export Expected A Posteriori person estimates
        person = estimate.2PL$person$EAP
      # Calculate estimate-to-true correlations for the IRT condition
        t.x = cor(person, x)
      # Calculate estimate-to-outcome correlations for the IRT condition
        t.y = cor(person, x1)
      # Calculate the rank order of participants when ability is estimated using IRT
        r.IRT = rank(person)
      # Find the 75th percentile
        s.class = person > summary(person)[5]
        m = matrix(NA, n, 1)
        for(i in 1:n){
        m[i,] = class[i] - s.class[i]
        }
      # Calculate false positive rate
        f.pos.irt = 100 * sum(m == -1)/n
      # Calculate false negative rate
        f.neg.irt = 100 * sum(m == 1)/n
    })}}
# Sim - Normal ---------------------------------------------------------------------
getwd()
setwd(paste("C:/Users/Administrator/Desktop/2PL Sim - Normal/Condition ", 1, sep = ""))
Optimal.N(n = 500, prop = .10, nitem = 20,
        ALPHA = 0.80, accept.info = 0.20,
        accept.rmsd = 0.05, File.Number = 1, nsim = 2000)
setwd(paste("C:/Users/Administrator/Desktop/2PL Sim - Normal/Condition ", 4, sep = ""))
Optimal.N(n = 250, prop = .10, nitem = 20,
        ALPHA = 0.80, accept.info = 0.20,
        accept.rmsd = 0.05, File.Number = 2, nsim = 2000)
setwd(paste("C:/Users/Administrator/Desktop/2PL Sim - Normal/Condition ", 7, sep = ""))
Optimal.N(n = 100, prop = .10, nitem = 20,
        ALPHA = 0.80, accept.info = 0.20,
        accept.rmsd = 0.05, File.Number = 3, nsim = 2000)

setwd(paste("C:/Users/Administrator/Desktop/2PL Sim - Normal/Condition ", 2, sep = ""))
Optimal.N(n = 500, prop = .30, nitem = 20,
        ALPHA = 0.80, accept.info = 0.20,
        accept.rmsd = 0.05, File.Number = 1, nsim = 2000)
setwd(paste("C:/Users/Administrator/Desktop/2PL Sim - Normal/Condition ", 5, sep = ""))
Optimal.N(n = 250, prop = .30, nitem = 20,
        ALPHA = 0.80, accept.info = 0.20,
        accept.rmsd = 0.05, File.Number = 2, nsim = 2000)
setwd(paste("C:/Users/Administrator/Desktop/2PL Sim - Normal/Condition ", 8, sep = ""))
Optimal.N(n = 100, prop = .30, nitem = 20,
        ALPHA = 0.80, accept.info = 0.20,
        accept.rmsd = 0.05, File.Number = 3, nsim = 2000)

setwd(paste("C:/Users/Administrator/Desktop/2PL Sim - Normal/Condition ", 3, sep = ""))
Optimal.N(n = 500, prop = .50, nitem = 20,
        ALPHA = 0.80, accept.info = 0.20,
        accept.rmsd = 0.05, File.Number = 1, nsim = 2000)
setwd(paste("C:/Users/Administrator/Desktop/2PL Sim - Normal/Condition ", 6, sep = ""))
Optimal.N(n = 250, prop = .50, nitem = 20,
        ALPHA = 0.80, accept.info = 0.20,
        accept.rmsd = 0.05, File.Number = 2, nsim = 2000)
setwd(paste("C:/Users/Administrator/Desktop/2PL Sim - Normal/Condition ", 9, sep = ""))
Optimal.N(n = 100, prop = .50, nitem = 20,
        ALPHA = 0.80, accept.info = 0.20,
        accept.rmsd = 0.05, File.Number = 3, nsim = 2000)
