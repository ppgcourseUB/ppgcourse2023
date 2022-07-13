#Modified by Víctor Noguerales, 14th December 2019

##This R script edits the .loci file on the basis of "per site theta within taxon" and "pairwise divergence between taxa" estimates.
#Also, it removes loci that have less than a specific number of taxa, and also those that are not variable.
#It needs the pegas R library (theta.s function)

#It can be also used to calculate the observed "per site theta within taxon" and "pairwise divergence between taxa" values
#in order to adjust the priors for some BPP analysis as AOO, for instance.
#For that purpose, change the number of the below parameters to a very high value. That way you would be sure you are estimating
#"per site theta within taxon" and "pairwise divergence between taxa" without excluding any loci.
#Thus, this allows you calculating the observed "theta" and "pairwise divergence" in your whole data set.
#You might want to do it before and after chopping your sequences...

#For example:
#max.theta.intra <- 0.999; #change number here
#max.div.inter <- 0.999;   #change number here

#Also, change this value according to the sequence length (before or after chopping)
#loci.seq.length <- 110;   #How long are the aligned sequences


library(pegas);
chopped.loci <- scan('DericorysSinOut_c85d5m18sh15_Chopped.loci', what = 'character', sep = '\n'); #Change file name
break.lines <- grep('//', chopped.loci);
cat('This .loci file contains a total of ', length(break.lines), 'loci'); #Let you know how many loci were in the original .loci file
index.lines <- c(0, break.lines);

loci.seq.length <- 110;   #How long are the aligned sequences
max.theta.intra <- 0.999; #Change number here
max.div.inter <- 0.999;   #Change number here
null.data <- rep(0, length(break.lines)*3); #For creating a data matrix
stat.mat <- matrix(null.data, nrow = length(break.lines), ncol = 3); #Data matrix for summary statistics
colnames(stat.mat) <- c('N_taxa', 'max_theta', 'max_p_dis');

for(i in 1:length(break.lines)){
	start.line <- index.lines[i] + 1;
	end.line <- break.lines[i] - 1;
	temp.alignment <- chopped.loci[start.line:end.line];

	s.pattern <- '>([[:alpha:]]+)[[:digit:]]+[[:space:]]+(.+)'; #Change this pattern according to your samples names
	sp.vec <- gsub(s.pattern, '\\1', temp.alignment);
	seq.vec <- gsub(s.pattern, '\\2', temp.alignment);
	uni.sp.vec <- unique(sp.vec);
	stat.mat[i,1] <- length(uni.sp.vec);

	theta.vec <- NULL; #Calculate theta first
	for(iter in 1:length(uni.sp.vec)){
		temp.sp <- uni.sp.vec[iter];
		temp.sample.index <- grep(temp.sp, sp.vec);
		if(length(temp.sample.index) > 1){ #At least 2 individual in a taxon
			var.sites <- NULL;
			for(ii in 1:(length(temp.sample.index)-1)){
				tt.seq1 <- unlist(strsplit(seq.vec[temp.sample.index[ii]], ''));
				tt.seq2 <- unlist(strsplit(seq.vec[temp.sample.index[ii + 1]], ''));
				temp.var.sites <- which(tt.seq1 != tt.seq2 & tt.seq1 != 'N' & tt.seq1 != '-' & tt.seq2 != 'N' & tt.seq2 != '-');
				var.sites <- c(var.sites, temp.var.sites);
			}
			n.segregating.sites <- length(unique(var.sites));
			temp.theta <- theta.s(n.segregating.sites, length(temp.sample.index))/loci.seq.length;
			theta.vec <- c(theta.vec, temp.theta);
		}else{  #Only one individual for a given taxon, just fill in 0
			theta.vec <- c(theta.vec, 0);
		}
	}
	max.theta <- max(theta.vec);
	stat.mat[i, 2] <- max.theta;

	p.dis.vec <- NULL; #Pairwise distance across samples from the same or different taxa
	for(iii in 1:(length(seq.vec)-1)){
		ttt.seq1 <- unlist(strsplit(seq.vec[iii], ''));
		for(itt in (iii+1):length(seq.vec)){
			ttt.seq2 <- unlist(strsplit(seq.vec[itt], ''));
			dis.sites <- which(ttt.seq1 != ttt.seq2 & ttt.seq1 != 'N' & ttt.seq2 != 'N' & ttt.seq1 != '-' & ttt.seq2 != '-');
			p.dis <- length(dis.sites)/loci.seq.length;
			p.dis.vec <- c(p.dis.vec, p.dis);
			#cat(iii, '\t', itt, '\n');
		}
	}
	max.p.dis <- max(p.dis.vec);
	stat.mat[i, 3] <- max.p.dis;
	cat(i, '\n');
}

#Note: The estimation of these parameters could take loooong if your dataset contains thousands of loci...

#Write the output and change file name to your preference. Save it for future use
#This .csv file can be used in the "1_count_n_varsites_locfile.R" script to display the "maximum theta per locus" and "Maximum pairwise divergence"
#Have a look at the "04_Prior setting in BPP4.0.R" if you wish to set the priors in BPP according to these empirical estimates

write.table(stat.mat, file = 'stat_chopped_loci.csv', sep = ',', quote = FALSE, row.names = FALSE);
