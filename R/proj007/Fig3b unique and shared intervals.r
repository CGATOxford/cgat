#Hs liver unique
# Figure Hs_liverunique_TSSgene_numbers
venn.diagram(
	x = list(
		HsLiver = 1:9200,
		TSSgene = 6943:26051
		),
	filename = "Hs_liverunique_TSSgene_numbers.tiff",
col = "transparent",
	fill = c("blue", "green"),
	alpha = 0.5,
	label.col = c("darkblue", "white", "darkgreen"),
	cex = 2.0,
	fontfamily = "arial",
	fontface = "bold",
	);

# Figure Hs_liverunique_TSSgene
venn.diagram(
	x = list(
		HsLiver = 1:9200,
		TSSgene = 6943:26051
		),
	filename = "Hs_liverunique_TSSgene.tiff",
col = "transparent",
	fill = c("blue", "green"),
	alpha = 0.5,
	cex = 0.0,
	cat.cex = 0.0,
);

#Hs Testes unique
# Figure Hs_testesunique_TSSgene_numbers
venn.diagram(
	x = list(
		Hstestes = 1:15118,
		TSSgene = 12808:31907
		),
	filename = "Hs_testesunique_TSSgene_numbers.tiff",
col = "transparent",
	fill = c("blue", "red"),
	alpha = 0.5,
	label.col = c("darkblue", "white", "darkgreen"),
	cex = 2.0,
	fontfamily = "arial",
	fontface = "bold",
	);

# Figure Hs_testesunique_TSSgene
venn.diagram(
	x = list(
		Hstestes = 1:15118,
		TSSgene = 12808:31907
		),
	filename = "Hs_testesunique_TSSgene.tiff",
col = "transparent",
	fill = c("blue", "red"),
	alpha = 0.5,
	cex = 0.0,
	cat.cex = 0.0,
);

#Hs Testes and Liver shared
# Figure Hs_testeslivershared_TSSgene_numbers
venn.diagram(
	x = list(
		Hsshared = 1:20847,
		TSSgene = 8031:26277
		),
	filename = "Hs_testeslivershared_TSSgene_numbers.tiff",
col = "transparent",
	fill = c("blue", "orange"),
	alpha = 0.5,
	label.col = c("darkblue", "white", "darkgreen"),
	cex = 2.0,
	fontfamily = "arial",
	fontface = "bold",
	);

# Figure Hs_testeslivershared_TSSgene
venn.diagram(
	x = list(
		Hsshared = 1:20847,
		TSSgene = 8031:26277
		),
	filename = "Hs_testeslivershared_TSSgene.tiff",
col = "transparent",
	fill = c("blue", "orange"),
	alpha = 0.5,
	cex = 0.0,
	cat.cex = 0.0,
);

###################################

