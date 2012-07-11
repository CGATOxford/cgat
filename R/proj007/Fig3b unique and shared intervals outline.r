#Hs liver unique
# Figure Hs_liverunique_TSSgene_numbers_line
venn.diagram(
	x = list(
		HsLiver = 1:9200,
		TSSgene = 6943:26051
		),
	filename = "Hs_liverunique_TSSgene_numbers_line.tiff",
col = "#58595B",
	fill = c("transparent"),
	alpha = 0.5,
	label.col = c("darkblue", "black", "darkgreen"),
	cex = 2.0,
	fontfamily = "arial",
	fontface = "bold",
	);

# Figure Hs_liverunique_TSSgene_line
venn.diagram(
	x = list(
		HsLiver = 1:9200,
		TSSgene = 6943:26051
		),
	filename = "Hs_liverunique_TSSgene_line.tiff",
col = "#58595B",
	fill = c("transparent"),
	alpha = 0.75,
	cex = 0.0,
	cat.cex = 0.0,
);

#Hs Testes unique
# Figure Hs_testesunique_TSSgene_numbers_line
venn.diagram(
	x = list(
		Hstestes = 1:15118,
		TSSgene = 12808:31907
		),
	filename = "Hs_testesunique_TSSgene_numbers_line.tiff",
col = "#58595B",
	fill = c("transparent"),
	alpha = 0.5,
	label.col = c("darkblue", "black", "darkgreen"),
	cex = 2.0,
	fontfamily = "arial",
	fontface = "bold",
	);

# Figure Hs_testesunique_TSSgene_line
venn.diagram(
	x = list(
		Hstestes = 1:15118,
		TSSgene = 12808:31907
		),
	filename = "Hs_testesunique_TSSgene_line.tiff",
col = "#58595B",
	fill = c("transparent"),
	alpha = 0.75,
	cex = 0.0,
	cat.cex = 0.0,
);

#Hs Testes and Liver shared
# Figure Hs_testeslivershared_TSSgene_numbers_line
venn.diagram(
	x = list(
		Hsshared = 1:20847,
		TSSgene = 8031:26277
		),
	filename = "Hs_testeslivershared_TSSgene_numbers_line.tiff",
col = "#58595B",
	fill = c("transparent"),
	alpha = 0.5,
	label.col = c("darkblue", "black", "darkgreen"),
	cex = 2.0,
	fontfamily = "arial",
	fontface = "bold",
	);

# Figure Hs_testeslivershared_TSSgene_line
venn.diagram(
	x = list(
		Hsshared = 1:20847,
		TSSgene = 8031:26277
		),
	filename = "Hs_testeslivershared_TSSgene_line.tiff",
col = "#58595B",
	fill = c("transparent"),
	alpha = 0.75,
	cex = 0.0,
	cat.cex = 0.0,
);

###################################

