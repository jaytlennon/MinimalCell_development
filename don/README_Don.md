# MinimalCell


*Reviewer 5 was brought in to evaluate statistical issues. Here is their primary concern:*

"Regarding the validity of the findings, despite the best efforts of the authors, 
the data profiles don't really fit the models that they are using. For example, 
see Fig. S1, where the authors use the 4-parameter modified (Zwietering) Gompertz 
equation to estimate growth parameters. Although the model is suitable, one needs to 
be very careful in extracting summary statistics/parameters for different populations. 
It would help to show the variation in individual growth curves, the residuals on fitting 
them to the model, and the actual distributions of the underlying parameters. 
The fit on the yield is fine, but the lag time is an issue (with marginal R^2 of 0.65 in the 
dangerously flexible GLMM). I am concerned about the existence of outliers and the lack of 
more information about the averaging, the actual growth measurements for each replicate, 
among others. I find that information is missing, or at least it is not clear; for 
example in figure S1, what are the blue and red dots? What are the units for yield 
and lag time? If they are CV, why are they not in %, as it was for the max growth rate? 
Why we can't see the error bars, perhaps on a separate plot, and/or 
overlapping growth curves?)."

So first thing, R5 is confused about interpretation of Fig. S1, which is not a fit of the 
Gompertz. Here is how I've been thinking we can respond:


1. Point them to Figshare documents containing fits for all ~120 growth curves

2. Add to this a diagnostic analysis of the growth curves

	a. Maybe add a four-panel figure like this: 
	https://www.r-bloggers.com/2019/11/the-hidden-diagnostic-plots-for-the-lm-object/
	Problem is that this is not as easy with mle2(), but pasted code below as a start
	
	b. Extract parameters associated with diagnostics commonly applied to these models,
	including residual vs. predicted, normal quantile plots, autocorrelation of residuals, 
	and outliers. I have pasted code below. 
	
	c. With the diagnostic summary stats, make distributions showing in fact that the 
	Gompertz does a very good job of fitting the data. 
	
3. More rigorously analyze the GLMM that uses the Gompertz paramaters (i.e., umax, lag, 
and yield). 

	a. There are some outliers here and I think R5 is giving license to remove altogether.
	
	b. Run diagnostics (above), but not sure exactly how to do in context of mixed model.
	
	c. Statistically interpret contribution of random effect to overal model
	
	
***** Example code for diagnostics statistics to include in Gompertz code*****	

	a. Slope of residuals vs. fitted should not be different from zero. Test this with 
	p value of lm resids ~ fitted values; should be >0.05 
	# results$eResid[i] <-summary(lm(residuals(best.f1)~predict(best.f1)))$coefficients[8]
	
	b. Normal quantile plot for test of normality. Slope of standardized residuals to 
	theoretical quantiles should be 1.0.
	# results$norm[i] <-summary(lm(residuals(best.f1)~ (qnorm(ppoints(length(residuals(best.f1))), 
	mean(residuals(best.f1)), sd(residuals(best.f1))))))$r.squared 
	
	c.Scale-location. Use Durbin-Watson test to look for autocorrleaiton:
	# results$DW[i] <- dwtest(residuals(best.f1) ~ 1)$p.value
	
	d. Residuals vs. leverage. More or less an attempt to identify influential i.e., 
	outlier residuals, which we can test with z-scores
	# results$z[i] <-sum(abs((residuals(fit4) - mean(residuals(fit4))) / sd(residuals(fit4)))>3)
	
	e. Report is "goodness of fit" with root mean square error:
	# results$RSME[i] <- sqrt(mean(residuals(best.f1))^2)
	


***** Code to make mle2 diagnostic four-panel figure*****	


# par(mfrow = c(2, 2), mar = c(5.1, 4.1, 4.1, 2.1))
#
# panel 1 = Residuals vs. fitted
# plot(predict(fit2), residuals(fit2), xlab = "Fitted values", ylab = "Residuals", las = 1)
# abline(h = 0, lty = 2, lwd = 1.5, col = "red")
#
# panel 2 = Normal Q-Q fit
# t_quant <- qnorm(ppoints(length(residuals(fit2))), mean(residuals(fit2)), sd(residuals(fit2)))
# plot(t_quant, residuals(fit2), main = "Normal Q-Q Plot", xlab = "Theoretical Quantiles", ylab = "Residuals")
# abline(a = mean(residuals(fit2)), b = sd(residuals(fit2)), col = "red", lwd = 1.5, lty = 2)
# 
# panel 3 = Scale location
# stdres <- residuals(fit2)/sqrt(coef(fit2)[5]) # need to specify error for last term
# sqrt_stdres <- sqrt(abs(stdres))
# plot(predict(fit2), sqrt_stdres, main = "Scale-Location Plot", xlab = "Fitted values",  
#      ylab = "Square root of standardized residuals")
# abline(h = c(1, -1)*qnorm(0.975)/sqrt(nrow(resp)), lty = 2, col = "red")
#
# panel 4 = Residuals vs. leverage
# residuals <- resid(fit2)
# std_resid <- residuals/sd(residuals)
# fmla <- formula(fit2)