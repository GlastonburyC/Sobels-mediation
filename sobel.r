
# Requires:
# covariate file
# dosage for cis-SNP # as.matrix(as.numeric(cis_genotype))
# expression residuals for trans genes and cis gene


#Calculate the percent mediation according to how the beta changes on conditions for the cis transcript
percent_mediation<- function(beta_orig , beta_adj){

	mediation_score = (beta_orig - beta_adj)/beta_orig
	return(mediation_score)
}

# change this function depending on covariates and whether you're modelling interactions or main effects.
mediation.test = function (mv, iv, dv, age, age_sq, bmi) 
{
    if (any(is.na(mv))) 
        stop("Mediator contains missing value(s)")
    if (any(is.na(iv))) 
        stop("Mediator contains missing value(s)")
    if (any(is.na(dv))) 
        stop("Mediator contains missing value(s)")
    nm = length(mv)
    ni = length(iv)
    nd = length(dv)
    if (nm != ni | nm != nd | ni != nd) 
        stop("Variables have different lengths.")
    tmp = summary(lm(mv ~ iv*bmi + age + age_sq ))
    a = tmp$coef[2, 1]
    sa = tmp$coef[2, 2]
    tmp = summary(lm(dv ~ mv + iv*bmi + age + age_sq))
    b = tmp$coef[3, 1]
    sb = tmp$coef[3, 2]
    tmp1 = b^2 * sa^2 + a^2 * sb^2
    tmp2 = sa^2 * sb^2
    zsob = a * b/sqrt(tmp1)
    psob = pnorm(-abs(zsob)) * 2
    zaro = a * b/sqrt(tmp1 + tmp2)
    paro = pnorm(-abs(zaro)) * 2
    if (tmp1 > tmp2) {
        zgm = a * b/sqrt(tmp1 - tmp2)
        pgm = pnorm(-abs(zgm)) * 2
    }
    else {
        zgm = NA
        pgm = NA
    }
    p.value = c(psob, paro, pgm)
    z.value = c(zsob, zaro, zgm)
    out = data.frame(rbind(z.value, p.value))
    names(out) = c("Sobel", "Aroian", "Goodman")
    out
}

covs<-read.table("/Users/Craig/Copy/interaction_plots/covs_matrix.txt")
age<-as.matrix(as.numeric(covs[1,]))
age_sq<-as.matrix(as.numeric(covs[2,]))
bmi<-as.matrix(as.numeric(covs[3,]))



results<-NULL
results=data.frame()

for(i in 1:dim(trans_expression)[1]) { 
	trans_gene_expr = as.matrix(as.double(trans_expression[i,2:length(trans_expression)]))
	
	gene_name = trans_expression[i,1]
	exon = trans_expression[i,2]

	# for each trans gene, calculate how it's effect changes when conditioned on the local cis transcript.
	orig_model = lm(trans_gene_expr ~ age + age_sq + bmi + cis_genotype*bmi)
	adj_model = lm (trans_gene_expr ~ age + age_sq + bmi + cis_expression + cis_genotype*bmi)

	orig_pval=summary(orig_model)$coefficients[6,4]
	adj_pval=summary(adj_model)$coefficients[7,4]
	# Extract the original and conditioned beta for the interaction
	beta_orig<-summary(orig_model)$coefficients[6,1]
	beta_adj<-summary(adj_model)$coefficients[7,1]

	# Call the function to calculate cis-mediation proportion
	med_score=percent_mediation(beta_orig,beta_adj)

	# Estimate p-value using sobel's test of mediation
	p_value= mediation.test(cis_expression,cis_genotype,trans_gene_expr,age,age_sq,bmi)[2,1]

	results = rbind(results, data.frame(gene_name ,exon, med_score , p_value, beta_orig, beta_adj,orig_pval,adj_pval))
}


pdf("ALG9_cis_mediation_analysis_boncorrected.pdf")
 p<-ggplot(results,aes(med_score,-log10(p_value)))
 p + geom_point() + geom_vline(h=0) + ggtitle("ALG9 trans network: Sobel's test of mediation") + theme_bw() + geom_hline(yintercept=2.97) + xlab("Mediation proportion (B-Badj/B)")
 m<-ggplot(results,aes(-log10(orig_pval),-log10(adj_pval)))
 m + geom_point() + xlab("ALG9 SNP association -log10(p-value)") + ylab("conditioned ALG9 SNP association -log10(p-value)") + ggtitle("conditioned on ALG9 expression") + theme_bw()
 s<-ggplot(results,aes(beta_adj,beta_orig))
 s + geom_point() + xlab("Conditioned trans gene beta for ALG9") + ylab("Original trans gene beta") + ggtitle("conditioned on ALG9 expression") + theme_bw()
dev.off()



