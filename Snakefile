mes_chr = list(range(1,23))+["X"]
HIC_DIR  = "/home/rochevin/Documents/PROJET_INGE/HiC_BENJ/data/"
HIC_FILES={i:HIC_DIR+i+".hic" for i in ["siCTRL_DIVA","siCTRL_OHT","siPER2_DIVA","siPER2_OHT"]}
#mes_res = "2500000 1000000 500000 10000".split(" ")
mes_res = "1000000 500000 250000 100000 50000 25000 10000 5000 1000".split(" ")
mes_res = {str(round(int(i)/1000))+"kb":int(i) for i in mes_res}
JUICER = "/mnt/NAS/DATA/HIGH_THROUGHPUT_GENOMICS_DIvA/Hi-C/HiC_F_Ben/juicer/CPU/common/juicer_tools.jar"
HTC_FORMAT = "create_HiTC_snakeformat.R"

def getFILES(wildcards):
	return(HIC_FILES[wildcards.exp])
def getRES(wildcards):
	return(mes_res[wildcards.reskb])

rule all:
	input:
		# expand("data/{norm}/dump_{obserOE}_{norm}_{exp}_chr{chr1}_chr{chr2}_{reskb}.txt.gz",
		# 	norm = ["KR"],
		# 	exp=HIC_FILES.keys(),
		# 	chr1 = mes_chr,chr2 = mes_chr,
		# 	reskb = mes_res.keys()),
		expand("TRANS/{norm}/HiTC/{obserOE}/HTC_{obserOE}_{norm}_{exp}_full_matrix_{reskb}.rds",
			norm = ["KR"],
			obserOE=["observed"],
			exp=HIC_FILES.keys(),
			reskb = mes_res.keys()),
		expand("TRANS/{norm}/dump_{obserOE}_{norm}_{exp}_chr{chr1}_chr{chr2}_{reskb}.txt.gz",
			norm = ["KR"],
			obserOE=["oe","observed"],
			exp=HIC_FILES.keys(),
			reskb = ["100kb","500kb","1000kb"],chr1 = mes_chr,chr2 = mes_chr)


rule dump_matrix:
	input:
		getFILES
	output:
		"TRANS/{norm}/dump_{obserOE}_{norm}_{exp}_chr{chr1}_chr{chr2}_{reskb}.txt.gz"
	params:
		my_res=getRES,
		juicer=JUICER,
		norm = "{norm}",
		chr1 = "{chr1}",
		chr2 = "{chr2}",
		obserOE="{obserOE}",
		out="TRANS/{norm}/dump_{obserOE}_{norm}_{exp}_chr{chr1}_chr{chr2}_{reskb}.txt"
	shell:
		"java -jar {params.juicer} "
		"dump {params.obserOE} "
		"{params.norm} {input} "
		"chr{params.chr1} chr{params.chr2} "
		"BP {params.my_res} "
		"{params.out} && gzip {params.out}"


rule HiTC_format:
	input:
		expand("TRANS/{{norm}}/dump_{{obserOE}}_{{norm}}_{{exp}}_chr{chr1}_chr{chr2}_{{reskb}}.txt.gz",chr1 = mes_chr,chr2 = mes_chr)
	output:
		"TRANS/{norm}/HiTC/{obserOE}/HTC_{obserOE}_{norm}_{exp}_full_matrix_{reskb}.rds"
	params:
		dumpdir = "TRANS/{norm}/",
		norm = "{norm}",
		exp = "{exp}",
		obserOE="{obserOE}",
		my_res=getRES,
		whichorder="custom"
	script:
		HTC_FORMAT

