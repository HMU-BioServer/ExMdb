<%@ page language="java" import="java.io.*,java.sql.*,java.util.*,java.text.DecimalFormat,wp.base.*" contentType="text/html; charset=UTF-8" pageEncoding="UTF-8"%>
<!doctype html>
<html lang="zxx">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<link rel="stylesheet" href="assets/css/bootstrap.min.css">
<link rel="stylesheet" href="assets/css/animate.min.css">
<link rel="stylesheet" href="assets/fonts/flaticon.css">
<link rel="stylesheet" href="assets/css/boxicons.min.css">
<link rel="stylesheet" href="assets/css/owl.carousel.min.css">
<link rel="stylesheet" href="assets/css/owl.theme.default.min.css">
<link rel="stylesheet" href="assets/css/magnific-popup.css">
<link rel="stylesheet" href="assets/css/nice-select.min.css">
<link rel="stylesheet" href="assets/css/meanmenu.css">
<link rel="stylesheet" href="assets/css/style.css">
<link rel="stylesheet" href="assets/css/responsive.css">
<link rel="icon" type="image/png" href="assets/images/favicon.png">

<link rel="stylesheet" type="text/css" href="select/demo.css"/>
<link rel="stylesheet" type="text/css" href="select/style-adsila.css" />
<link rel="stylesheet" href="select/selectpage_cell.css" type="text/css">

<title>ExMdb : Pseudotime pathway</title>



<style>
.search_item{
	font-size:20px;
}

.imm_select{
	width:100%;
}


.theme-btn {
	margin-top:70%;
	float:left;
	padding:10%;
	font-size: 15px;
}

.graph_qi {
    position: fixed;
    top: 0;
    left: 0;
    width: 100%;
    height: 100%;
    z-index: 99999;
    background: #00000073;
    display:none;
}

div.sp_clear_btn {
    font-size: 25px;
}

ul.sp_element_box{
	padding: 5px 0 0 3px;
}

.sp_result_area{
	width:500px
}

.sp_results > li{
	font-size: 10px;
}

.list{
	height: 200px;
    overflow-y: auto !important;
}

.score{
	width:100%;
}

.score > ul{
	height:120px
}
</style>

	<style type="text/css">
		/* NOTE: The styles were added inline because Prefixfree needs access to your styles and they must be inlined if they are on local disk! */
		/* =========
		Get Fonts */
		/* ================
		   Assign Variables */
		/* ===========================
		   Setup Mixins/Helper Classes */
		.clearfix:after, .tab-nav-small:after {
		  content: ".";
		  display: block;
		  height: 0;
		  clear: both;
		  visibility: hidden;
		}

		/* ==========
		   Setup Page */


		/* ===========
		   Tab Styling */
		.tab-group {
		  position: relative;
		  border: 1px solid #eee;
		  margin-top: 2.5em;
		  border-radius: 0 0 10px 10px;
		}
		.tab-group section {
		  opacity: 0;
		  height: 0;
		  padding: 0 1em;
		  overflow: hidden;
		  transition: opacity 0.4s ease, height 0.4s ease;
		}
		.tab-group section.active {
		  opacity: 1;
		  height: auto;
		  overflow: visible;
		}

		.tab-nav-small {
		  list-style: none;
		  margin: -2.5em -1px 0 0;
		  padding: 0;
		  height: 2.5em;
		  overflow: hidden;
		}
		.tab-nav-small li {
		  display: inline;
		}
		.tab-nav-small li a {
		  top: 1px;
		  position: relative;
		  display: block;
		  float: left;
		  border-radius: 10px 10px 0 0;
		  background: #eee;
		  line-height: 2em;
		  padding: 0 1em;
		  text-decoration: none;
		  color: grey;
		  margin-top: .5em;
		  margin-right: 1px;
		  transition: background .2s ease, line-height .2s ease, margin .2s ease;
		}
		.tab-nav-small li.active a {
		  background: #6671c1;
		  color: white;
		  line-height: 2.5em;
		  margin-top: 0;
		}
	</style>

</head>
<body>


<div class="preloader">
  <div class="d-table">
    <div class="d-table-cell">
      <div class="spinner"></div>
    </div>
  </div>
</div>

<div class="graph_qi" name="graph_qi" id="graph_qi">
  <div class="d-table">
    <div class="d-table-cell">
      <div class="spinner"></div>
    </div>
  </div>
</div>


<div class="navbar-area" data-step="1" data-intro="This is a tooltip!">
  <div class="mobile-nav"><a href="index.html" class="logo"><img src="assets/images/logos/logo-1.png" alt="Logo"></a></div>
  <div class="main-nav">
    <div class="container">
      <nav class="navbar navbar-expand-md navbar-light "><a class="navbar-brand" href="index.html"><img src="assets/images/logos/logo-1.png" alt="Logo"></a>
        <div class="collapse navbar-collapse mean-menu" id="navbarSupportedContent">
          <ul class="navbar-nav m-auto" style="z-index:100">
            <li class="nav-item"><a href="index.html" class="nav-link">Home</a></li>
            <li class="nav-item"><a href="exmdb_browse.jsp" target="_blank" class="nav-link">Browse</a></li>
            <li class="nav-item"><a href="#" class="nav-link">Search <i class='bx bx-caret-down'></i></a>
              <ul class="dropdown-menu">
                <li class="nav-item"><a href="exmdb_search_highput.jsp" class="nav-link" target="_blank">Search [ High throughput ]</a></li>
                <li class="nav-item"><a href="exmdb_search_exp.jsp" class="nav-link" target="_blank">Search [ Experimental validation ]</a></li>
                <li class="nav-item"><a href="exmdb_search_biomarker.jsp" class="nav-link" target="_blank">Search [ Cancer Biomarkers ]</a></li>
                <li class="nav-item"><a href="exmdb_search_scdataset.jsp" class="nav-link" target="_blank">Search [ Single-cell ]</a></li>
              </ul>
            </li>
            <li class="nav-item"><a href="#" class="nav-link active">Tools <i class='bx bx-caret-down'></i></a>
              <ul class="dropdown-menu">
                <li class="nav-item"><a href="exmdb_cell_out.jsp" target="_blank" class="nav-link">ExMdb-Cellloaction </a></li>
                <li class="nav-item"><a href="exmdb_imm_out.jsp" target="_blank" class="nav-link">ExMdb-Immunity </a></li><li class="nav-item"><a href="#" class="nav-link">ExMdb-Survival Tools<i class='bx bx-caret-down'></i></a><ul class="dropdown-menu"><li class="nav-item"><a href="exmdb_survival_out.jsp" target="_blank" class="nav-link" style="text-transform: none;">mRNA-Survival Analyze</a></li><li class="nav-item"><a href="exmdb_survival_miRNA_out.jsp" target="_blank" class="nav-link" style="text-transform: none;">miRNA-Survival Analyze</a></li> <li class="nav-item"><a href="exmdb_survival_out.jsp" target="_blank" class="nav-link" style="text-transform: none;">lncRNA-Survival Analyze</a></li></ul></li>
                <li class="nav-item"><a href="exmdb_net_out.jsp" target="_blank" class="nav-link">ExMdb-Network </a></li><li class="nav-item"><a href="exmdb_blast_out.jsp" target="_blank" class="nav-link">ExMdb-BLAST </a></li><li class="nav-item"><a href="exmdb_function_time.jsp" target="_blank" class="nav-link">ExMdb-Pseudotime pathway</a></li><li class="nav-item"><a href="exmdb_single_vis.jsp" target="_blank" class="nav-link">ExMdb-Single Cell Visualization</a></li>
				<li class="nav-item"><a href="#" class="nav-link">ExMdb-FunctionEnrichment <i class='bx bx-caret-down'></i></a><ul class="dropdown-menu"><li class="nav-item"><a href="exmdb_function_mRNA_out.jsp" target="_blank" class="nav-link" style="text-transform: none;">mRNA-FunctionEnrichment</a></li><li class="nav-item"><a href="exmdb_function_out_miRNA.jsp" target="_blank" class="nav-link" style="text-transform: none;">miRNA-FunctionEnrichment</a></li> <li class="nav-item"><a href="exmdb_function_out.jsp" target="_blank" class="nav-link" style="text-transform: none;">lncRNA-FunctionEnrichment</a></li></ul></li>
              </ul>
            </li>


            <li class="nav-item"><a href="exmdb_statistic.jsp" target="_blank" class="nav-link">Statistic </a></li><li class="nav-item"><a href="exmdb_help.html" target="_blank" class="nav-link">Help </a></li>
            
            <li class="nav-item"><a href="exmdb_download.html" target="_blank" class="nav-link">Download </a></li><li class="nav-item"><a href="exmdb_contact.html" target="_blank" class="nav-link">Contact Us</a></li>
          </ul>
          
        </div>
      </nav>
    </div>
  </div>
  
</div>

<div class="inner-banner">
  <div class="container">
    <div class="inner-title text-center">
      <h3>Pseudotime pathway</h3>
      <ul>
        <li><a href="index.html">Home</a></li>
        <li><i class='bx bx-chevrons-right'></i></li>
        <li>Pseudotime pathway</li>
      </ul>
    </div>
  </div>
  <div class="inner-shape"><img src="assets/images/shape/inner-shape.png" alt="Images"></div>
</div>
  <div class="container wow fadeInUp animated" style="max-width: 1440px;padding-top: 2.5rem;">

		<div class="contact-header">
			<h1>Pseudotime pathway</h1>
		</div>	
		<div class="contact-wrapper">
	    
		    <p class="sub_title_function" style="margin-top:15px"><img src="img/iconfont/pse.png" style="margin-top: -7px;margin-right: 7px;">ExMdb-Pseudotime Pathway</p>
	    	<p class="sub_text_function">Using Pseudotime analysis to delineate cellular functional pathway changes at different time points.</p><br>
	    <div class="tab-group">
		    <section id="tab1" title="Biological process">
		    <form action="Code_test.jsp" target="exmdb_function"  method="post" role="form" id="imm_form" name="imm_form">
		    

		    
            <div class="row">
            <div class="col-md-1 col-sm-1 col-xs-1"></div>
            <div class="col-md-10 col-sm-10 col-xs-10">
							<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name"> Biological process</span>
									<span class="menu__item-label">(Enter pathway/keywords)</span>
								</a>
								<input type="text" name="searchname" id="searchname" class="form-control" value="GOBP_2_OXOGLUTARATE_METABOLIC_PROCESS,GOBP_3_PHOSPHOADENOSINE_5_PHOSPHOSULFATE_METABOLIC_PROCESS,GOBP_3_UTR_MEDIATED_MRNA_DESTABILIZATION,GOBP_3_UTR_MEDIATED_MRNA_STABILIZATION,GOBP_ACETYL_COA_BIOSYNTHETIC_PROCESS,GOBP_ACETYL_COA_BIOSYNTHETIC_PROCESS_FROM_PYRUVATE,GOBP_ACETYL_COA_METABOLIC_PROCESS,GOBP_ACID_SECRETION,GOBP_ACIDIC_AMINO_ACID_TRANSPORT" placeholder=" VERY LONG CHAIN FATTY ACID METABOLIC PROCESS / SMALL NUCLEOLAR RIBONUCLEOPROTEIN COMPLEX ASSEMBLY ..."  style="font-size:15px;padding:20px 15px">
							</nav>
            </div></div>
            <div class="row" style="padding-bottom: 2rem;">
            <div class="col-md-1 col-sm-1 col-xs-1"></div>
            <div class="col-md-5 col-sm-5 col-xs-5">
            				<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name">DATASET ID</span>
									<span class="menu__item-label">(GSE/COVID)</span>
								</a>
								<select class="language-list-item imm_select" name="gseid" id="gseid">
									<option value="BLCA_GSE145281_aPDL1">BLCA_GSE145281_aPDL1 [ N = 10 ] , Cell Count = [ 14474 ]</option>
									<option value="CLL_GSE111014">CLL_GSE111014 [ N = 4 ] , Cell Count = [ 30106 ]</option>
									<option value="KIRC_GSE145281_aPDL1">KIRC_GSE145281_aPDL1 [ N = 4 ] , Cell Count = [ 44220 ]</option>
									<option value="MCC_GSE117988_aPD1aCTLA4">MCC_GSE117988_aPD1aCTLA4 [ N = 1 ] , Cell Count = [ 10148 ]</option>
									<option value="MCC_GSE118056_aPDL1">MCC_GSE118056_aPDL1 [ N = 1 ] , Cell Count = [ 11024 ]</option>
									<option value="NSCLC_GSE127471">NSCLC_GSE127471 [ N = 1 ] , Cell Count = [ 1108 ]</option>
									<option value="PBMC_30K_10X">PBMC_30K_10X [ N = 1 ] , Cell Count = [ 29079 ]</option>
									<option value="PBMC_8K_10X">PBMC_8K_10X [ N = 1 ] , Cell Count = [ 8488 ]</option>
									<option value="PC_GSE67980">PC_GSE67980 [ N = 13 ] , Cell Count = [ 117 ]</option>
									<option value="MEL_GSE157743">MEL_GSE157743 [ N = 22 ] , Cell Count = [ 102 ]</option>
									<option value="COVID19_GSE168212">COVID19_GSE168212 [ N = 1 ] , Cell Count = [ 4935 ]</option>
									<option value="COVID19_GSE166992">COVID19_GSE166992 [ N = 9 ] , Cell Count = [ 20348 ]</option>
									<option value="COVID19_GSE154567_severe">COVID19_GSE154567_severe [ N = 3 ] , Cell Count = [ 15190 ]</option>
									<option value="COVID19_GSE154567_moderate">COVID19_GSE154567_moderate [ N = 3 ] , Cell Count = [ 22809 ]</option>
									<option value="COVID19_GSE154567_recovery">COVID19_GSE154567_recovery [ N = 3 ] , Cell Count = [ 29275 ]</option>
									<option value="COVID19_GSE167118_BALF_moderate">COVID19_GSE167118_BALF_moderate [ N = 2 ] , Cell Count = [ 3016 ]</option>
									<option value="COVID19_GSE167118_BALF_severe">COVID19_GSE167118_BALF_severe [ N = 6 ] , Cell Count = [ 23581 ]</option>
									<option value="COVID19_GSE167118_blood_moderate">COVID19_GSE167118_blood_moderate [ N = 1 ] , Cell Count = [ 7270 ]</option>
									<option value="COVID19_GSE167118_blood_severe">COVID19_GSE167118_blood_severe [ N = 6 ] , Cell Count = [ 54094 ]</option>
									<option value="BP_GSE167118_BALF_moderate">BP_GSE167118_BALF_moderate [ N = 3 ] , Cell Count = [ 5560 ]</option>
									<option value="BP_GSE167118_BALF_no">BP_GSE167118_BALF_no [ N = 1 ] , Cell Count = [ 3349 ]</option>
									<option value="BP_GSE167118_blood_moderate">BP_GSE167118_blood_moderate [ N = 3 ] , Cell Count = [ 14866 ]</option>
									<option value="BP_GSE167118_blood_no">BP_GSE167118_blood_no [ N = 1 ] , Cell Count = [ 6988 ]</option>
									<option value="COVID19_GSE169503">COVID19_GSE169503 [ N = 4 ] , Cell Count = [ 4196 ]</option>
									<option value="SLE_GSE142016">SLE_GSE142016 [ N = 3 ] , Cell Count = [ 18531 ]</option>
									<option value="LIHC_GSE107747">LIHC_GSE107747 [ N = 2 ] , Cell Count = [ 9950 ]</option>
									<option value="Lymphoma_GSE124899">Lymphoma_GSE124899 [ N = 2 ] , Cell Count = [ 15124 ]</option>
									<option value="IPEX_GSE167976">IPEX_GSE167976 [ N = 4 ] , Cell Count = [ 10375 ]</option>
					            </select>
							</nav>

            </div>
            
            <div class="col-md-2 col-sm-2 col-xs-2">
            				<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name">No. of Split</span>
									<span class="menu__item-label"></span>
								</a>
								<select class="language-list-item imm_select" name="top" id="top">
					              <option value="10" selected="selected">10</option>
					              <option value="20">20</option>
					              <option value="30">30</option>
					              <option value="40">40</option>
					              <option value="50">50</option>
					            </select>
							</nav>

            </div>
            <div class="col-md-3 col-sm-3 col-xs-3">
            				<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name">Score</span>
									<span class="menu__item-label"></span>
								</a>
								<select class="score" name="method" id="method">
					              <option value="normalize" selected="selected">Normalize [ X - min(X) ]</option>
					              <option value="raw">Raw Data</option>
					              <option value="zero">Zero Start [ Xn - X1 ]</option>
					            </select>
							</nav>
            </div>
            <div class="col-md-1 col-sm-1 col-xs-1">
           		<button class="default-btn btn-bg-two border-radius-50 theme-btn" id="submit_botton" type="submit">SUBMIT</button>
            </div>
            </div>
            </form>
		    </section>
		    <section id="tab2" title=" Cellular Component">
		    <form action="Code_test.jsp" target="exmdb_function"  method="post" role="form" id="imm_form" name="imm_form">
			<div class="row">
            <div class="col-md-1 col-sm-1 col-xs-1"></div>
            <div class="col-md-10 col-sm-10 col-xs-10">
							<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name"> Cellular Component</span>
									<span class="menu__item-label">(Enter pathway/keywords)</span>
								</a>
								<input type="text" name="searchname_cc" id="searchname_cc" class="form-control" value="GOCC_90S_PRERIBOSOME,GOCC_9PLUS0_NON_MOTILE_CILIUM,GOCC_9PLUS2_MOTILE_CILIUM,GOCC_A_BAND,GOCC_ACETYLCHOLINE_GATED_CHANNEL_COMPLEX,GOCC_ACROSOMAL_MEMBRANE,GOCC_ACROSOMAL_VESICLE,GOCC_ACTIN_BASED_CELL_PROJECTION,GOCC_ACTIN_CYTOSKELETON" placeholder=" NUCLEAR EXOSOME RNASE COMPLEX / VACUOLAR PROTON TRANSPORTING V TYPE ATPASE V0 DOMAIN ..."  style="font-size:15px;padding:20px 15px">
							</nav>
            </div></div>
            <div class="row" style="padding-bottom: 2rem;">
            
            
             <div class="col-md-1 col-sm-1 col-xs-1"></div>
            
            <div class="col-md-5 col-sm-5 col-xs-5">
            				<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name">DATASET ID</span>
									<span class="menu__item-label">(GSE/COVID)</span>
								</a>
								<select class="language-list-item imm_select" name="gseid" id="gseid">
									<option value="BLCA_GSE145281_aPDL1">BLCA_GSE145281_aPDL1 [ N = 10 ] , Cell Count = [ 14474 ]</option>
									<option value="CLL_GSE111014">CLL_GSE111014 [ N = 4 ] , Cell Count = [ 30106 ]</option>
									<option value="KIRC_GSE145281_aPDL1">KIRC_GSE145281_aPDL1 [ N = 4 ] , Cell Count = [ 44220 ]</option>
									<option value="MCC_GSE117988_aPD1aCTLA4">MCC_GSE117988_aPD1aCTLA4 [ N = 1 ] , Cell Count = [ 10148 ]</option>
									<option value="MCC_GSE118056_aPDL1">MCC_GSE118056_aPDL1 [ N = 1 ] , Cell Count = [ 11024 ]</option>
									<option value="NSCLC_GSE127471">NSCLC_GSE127471 [ N = 1 ] , Cell Count = [ 1108 ]</option>
									<option value="PBMC_30K_10X">PBMC_30K_10X [ N = 1 ] , Cell Count = [ 29079 ]</option>
									<option value="PBMC_8K_10X">PBMC_8K_10X [ N = 1 ] , Cell Count = [ 8488 ]</option>
									<option value="PC_GSE67980">PC_GSE67980 [ N = 13 ] , Cell Count = [ 117 ]</option>
									<option value="MEL_GSE157743">MEL_GSE157743 [ N = 22 ] , Cell Count = [ 102 ]</option>
									<option value="COVID19_GSE168212">COVID19_GSE168212 [ N = 1 ] , Cell Count = [ 4935 ]</option>
									<option value="COVID19_GSE166992">COVID19_GSE166992 [ N = 9 ] , Cell Count = [ 20348 ]</option>
									<option value="COVID19_GSE154567_severe">COVID19_GSE154567_severe [ N = 3 ] , Cell Count = [ 15190 ]</option>
									<option value="COVID19_GSE154567_moderate">COVID19_GSE154567_moderate [ N = 3 ] , Cell Count = [ 22809 ]</option>
									<option value="COVID19_GSE154567_recovery">COVID19_GSE154567_recovery [ N = 3 ] , Cell Count = [ 29275 ]</option>
									<option value="COVID19_GSE167118_BALF_moderate">COVID19_GSE167118_BALF_moderate [ N = 2 ] , Cell Count = [ 3016 ]</option>
									<option value="COVID19_GSE167118_BALF_severe">COVID19_GSE167118_BALF_severe [ N = 6 ] , Cell Count = [ 23581 ]</option>
									<option value="COVID19_GSE167118_blood_moderate">COVID19_GSE167118_blood_moderate [ N = 1 ] , Cell Count = [ 7270 ]</option>
									<option value="COVID19_GSE167118_blood_severe">COVID19_GSE167118_blood_severe [ N = 6 ] , Cell Count = [ 54094 ]</option>
									<option value="BP_GSE167118_BALF_moderate">BP_GSE167118_BALF_moderate [ N = 3 ] , Cell Count = [ 5560 ]</option>
									<option value="BP_GSE167118_BALF_no">BP_GSE167118_BALF_no [ N = 1 ] , Cell Count = [ 3349 ]</option>
									<option value="BP_GSE167118_blood_moderate">BP_GSE167118_blood_moderate [ N = 3 ] , Cell Count = [ 14866 ]</option>
									<option value="BP_GSE167118_blood_no">BP_GSE167118_blood_no [ N = 1 ] , Cell Count = [ 6988 ]</option>
									<option value="COVID19_GSE169503">COVID19_GSE169503 [ N = 4 ] , Cell Count = [ 4196 ]</option>
									<option value="SLE_GSE142016">SLE_GSE142016 [ N = 3 ] , Cell Count = [ 18531 ]</option>
									<option value="LIHC_GSE107747">LIHC_GSE107747 [ N = 2 ] , Cell Count = [ 9950 ]</option>
									<option value="Lymphoma_GSE124899">Lymphoma_GSE124899 [ N = 2 ] , Cell Count = [ 15124 ]</option>
									<option value="IPEX_GSE167976">IPEX_GSE167976 [ N = 4 ] , Cell Count = [ 10375 ]</option>
					            </select>
							</nav>

            </div>
            
            <div class="col-md-2 col-sm-2 col-xs-2">
            				<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name">No. of Split</span>
									<span class="menu__item-label"></span>
								</a>
								<select class="language-list-item imm_select" name="top" id="top">
					              <option value="10" selected="selected">10</option>
					              <option value="20">20</option>
					              <option value="30">30</option>
					              <option value="40">40</option>
					              <option value="50">50</option>
					            </select>
							</nav>

            </div>     
            <div class="col-md-3 col-sm-3 col-xs-3">
            				<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name">Score</span>
									<span class="menu__item-label"></span>
								</a>
								<select class="score" name="method" id="method">
					              <option value="normalize" selected="selected">Normalize [ X - min(X) ]</option>
					              <option value="raw">Raw Data</option>
					              <option value="zero">Zero Start [ Xn - X1 ]</option>
					            </select>
							</nav>
            </div>
            <div class="col-md-1 col-sm-1 col-xs-1">
           		<button class="default-btn btn-bg-two border-radius-50 theme-btn" id="submit_botton" type="submit">SUBMIT</button>
            </div>
            </div>
            </form>
		    </section>
		    <section id="tab3" title="Molecular Function">
		    <form action="Code_test.jsp" target="exmdb_function"  method="post" role="form" id="imm_form" name="imm_form">
			<div class="row">
            <div class="col-md-1 col-sm-1 col-xs-1"></div>
            <div class="col-md-10 col-sm-10 col-xs-10">
							<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name"> Molecular Function</span>
									<span class="menu__item-label">(Enter pathway/keywords)</span>
								</a>
								<input type="text" name="searchname_mf" id="searchname_mf" class="form-control" value="GOMF_1_PHOSPHATIDYLINOSITOL_3_KINASE_ACTIVITY,GOMF_1_PHOSPHATIDYLINOSITOL_3_KINASE_REGULATOR_ACTIVITY,GOMF_1_PHOSPHATIDYLINOSITOL_BINDING,GOMF_14_3_3_PROTEIN_BINDING,GOMF_17_BETA_HYDROXYSTEROID_DEHYDROGENASE_NADPLUS_ACTIVITY,GOMF_17_BETA_HYDROXYSTEROID_DEHYDROGENASE_NADPPLUS_ACTIVITY,GOMF_2_IRON_2_SULFUR_CLUSTER_BINDING,GOMF_2_OXOGLUTARATE_DEPENDENT_DIOXYGENASE_ACTIVITY,GOMF_3_5_CYCLIC_AMP_PHOSPHODIESTERASE_ACTIVITY" placeholder=" ADENYL_NUCLEOTIDE_EXCHANGE_FACTOR_ACTIVITY / DNA_BINDING_TRANSCRIPTION_ACTIVATOR_ACTIVITY ..."  style="font-size:15px;padding:20px 15px">
							</nav>
            </div></div>
            <div class="row" style="padding-bottom: 2rem;">
            <div class="col-md-1 col-sm-1 col-xs-1"></div>
            <div class="col-md-5 col-sm-5 col-xs-5">
            				<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name">DATASET ID</span>
									<span class="menu__item-label">(GSE/COVID)</span>
								</a>
								<select class="language-list-item imm_select" name="gseid" id="gseid">
									<option value="BLCA_GSE145281_aPDL1">BLCA_GSE145281_aPDL1 [ N = 10 ] , Cell Count = [ 14474 ]</option>
									<option value="CLL_GSE111014">CLL_GSE111014 [ N = 4 ] , Cell Count = [ 30106 ]</option>
									<option value="KIRC_GSE145281_aPDL1">KIRC_GSE145281_aPDL1 [ N = 4 ] , Cell Count = [ 44220 ]</option>
									<option value="MCC_GSE117988_aPD1aCTLA4">MCC_GSE117988_aPD1aCTLA4 [ N = 1 ] , Cell Count = [ 10148 ]</option>
									<option value="MCC_GSE118056_aPDL1">MCC_GSE118056_aPDL1 [ N = 1 ] , Cell Count = [ 11024 ]</option>
									<option value="NSCLC_GSE127471">NSCLC_GSE127471 [ N = 1 ] , Cell Count = [ 1108 ]</option>
									<option value="PBMC_30K_10X">PBMC_30K_10X [ N = 1 ] , Cell Count = [ 29079 ]</option>
									<option value="PBMC_8K_10X">PBMC_8K_10X [ N = 1 ] , Cell Count = [ 8488 ]</option>
									<option value="PC_GSE67980">PC_GSE67980 [ N = 13 ] , Cell Count = [ 117 ]</option>
									<option value="MEL_GSE157743">MEL_GSE157743 [ N = 22 ] , Cell Count = [ 102 ]</option>
									<option value="COVID19_GSE168212">COVID19_GSE168212 [ N = 1 ] , Cell Count = [ 4935 ]</option>
									<option value="COVID19_GSE166992">COVID19_GSE166992 [ N = 9 ] , Cell Count = [ 20348 ]</option>
									<option value="COVID19_GSE154567_severe">COVID19_GSE154567_severe [ N = 3 ] , Cell Count = [ 15190 ]</option>
									<option value="COVID19_GSE154567_moderate">COVID19_GSE154567_moderate [ N = 3 ] , Cell Count = [ 22809 ]</option>
									<option value="COVID19_GSE154567_recovery">COVID19_GSE154567_recovery [ N = 3 ] , Cell Count = [ 29275 ]</option>
									<option value="COVID19_GSE167118_BALF_moderate">COVID19_GSE167118_BALF_moderate [ N = 2 ] , Cell Count = [ 3016 ]</option>
									<option value="COVID19_GSE167118_BALF_severe">COVID19_GSE167118_BALF_severe [ N = 6 ] , Cell Count = [ 23581 ]</option>
									<option value="COVID19_GSE167118_blood_moderate">COVID19_GSE167118_blood_moderate [ N = 1 ] , Cell Count = [ 7270 ]</option>
									<option value="COVID19_GSE167118_blood_severe">COVID19_GSE167118_blood_severe [ N = 6 ] , Cell Count = [ 54094 ]</option>
									<option value="BP_GSE167118_BALF_moderate">BP_GSE167118_BALF_moderate [ N = 3 ] , Cell Count = [ 5560 ]</option>
									<option value="BP_GSE167118_BALF_no">BP_GSE167118_BALF_no [ N = 1 ] , Cell Count = [ 3349 ]</option>
									<option value="BP_GSE167118_blood_moderate">BP_GSE167118_blood_moderate [ N = 3 ] , Cell Count = [ 14866 ]</option>
									<option value="BP_GSE167118_blood_no">BP_GSE167118_blood_no [ N = 1 ] , Cell Count = [ 6988 ]</option>
									<option value="COVID19_GSE169503">COVID19_GSE169503 [ N = 4 ] , Cell Count = [ 4196 ]</option>
									<option value="SLE_GSE142016">SLE_GSE142016 [ N = 3 ] , Cell Count = [ 18531 ]</option>
									<option value="LIHC_GSE107747">LIHC_GSE107747 [ N = 2 ] , Cell Count = [ 9950 ]</option>
									<option value="Lymphoma_GSE124899">Lymphoma_GSE124899 [ N = 2 ] , Cell Count = [ 15124 ]</option>
									<option value="IPEX_GSE167976">IPEX_GSE167976 [ N = 4 ] , Cell Count = [ 10375 ]</option>
					            </select>
							</nav>

            </div>
            
            <div class="col-md-2 col-sm-2 col-xs-2">
            				<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name">No. of Split</span>
									<span class="menu__item-label"></span>
								</a>
								<select class="language-list-item imm_select" name="top" id="top">
					              <option value="10" selected="selected">10</option>
					              <option value="20">20</option>
					              <option value="30">30</option>
					              <option value="40">40</option>
					              <option value="50">50</option>
					            </select>
							</nav>

            </div>     
            <div class="col-md-3 col-sm-3 col-xs-3">
            				<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name">Score</span>
									<span class="menu__item-label"></span>
								</a>
								<select class="score" name="method" id="method">
					              <option value="normalize" selected="selected">Normalize [ X - min(X) ]</option>
					              <option value="raw">Raw Data</option>
					              <option value="zero">Zero Start [ Xn - X1 ]</option>
					            </select>
							</nav>
            </div>
            <div class="col-md-1 col-sm-1 col-xs-1">
           		<button class="default-btn btn-bg-two border-radius-50 theme-btn" id="submit_botton" type="submit">SUBMIT</button>     	
            </div>
            </div></form>
		    </section>
		    <section id="tab4" title="50 Hallmaker">
		    <form action="Code_test.jsp" target="exmdb_function"  method="post" role="form" id="imm_form" name="imm_form">
			<div class="row">
            <div class="col-md-1 col-sm-1 col-xs-1"></div>
            <div class="col-md-10 col-sm-10 col-xs-10">
							<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name"> 50 Hallmaker</span>
									<span class="menu__item-label">(Enter pathway/keywords)</span>
								</a>
								<input type="text" name="searchname_hallmaker" id="searchname_hallmaker" class="form-control" value="HALLMARK_ADIPOGENESIS,HALLMARK_ALLOGRAFT_REJECTION,HALLMARK_ANDROGEN_RESPONSE,HALLMARK_ANGIOGENESIS,HALLMARK_APICAL_JUNCTION,HALLMARK_APICAL_SURFACE,HALLMARK_APOPTOSIS,HALLMARK_BILE_ACID_METABOLISM,HALLMARK_CHOLESTEROL_HOMEOSTASIS" placeholder=" INTERFERON_ALPHA_RESPONSE / ANGIOGENESIS ..."  style="font-size:15px;padding:20px 15px">
							</nav>
            </div></div>
            <div class="row" style="padding-bottom: 2rem;">
            <div class="col-md-1 col-sm-1 col-xs-1"></div>
            <div class="col-md-5 col-sm-5 col-xs-5">
            				<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name">DATASET ID</span>
									<span class="menu__item-label">(GSE/COVID)</span>
								</a>
								<select class="language-list-item imm_select" name="gseid" id="gseid">
									<option value="BLCA_GSE145281_aPDL1">BLCA_GSE145281_aPDL1 [ N = 10 ] , Cell Count = [ 14474 ]</option>
									<option value="CLL_GSE111014">CLL_GSE111014 [ N = 4 ] , Cell Count = [ 30106 ]</option>
									<option value="KIRC_GSE145281_aPDL1">KIRC_GSE145281_aPDL1 [ N = 4 ] , Cell Count = [ 44220 ]</option>
									<option value="MCC_GSE117988_aPD1aCTLA4">MCC_GSE117988_aPD1aCTLA4 [ N = 1 ] , Cell Count = [ 10148 ]</option>
									<option value="MCC_GSE118056_aPDL1">MCC_GSE118056_aPDL1 [ N = 1 ] , Cell Count = [ 11024 ]</option>
									<option value="NSCLC_GSE127471">NSCLC_GSE127471 [ N = 1 ] , Cell Count = [ 1108 ]</option>
									<option value="PBMC_30K_10X">PBMC_30K_10X [ N = 1 ] , Cell Count = [ 29079 ]</option>
									<option value="PBMC_8K_10X">PBMC_8K_10X [ N = 1 ] , Cell Count = [ 8488 ]</option>
									<option value="PC_GSE67980">PC_GSE67980 [ N = 13 ] , Cell Count = [ 117 ]</option>
									<option value="MEL_GSE157743">MEL_GSE157743 [ N = 22 ] , Cell Count = [ 102 ]</option>
									<option value="COVID19_GSE168212">COVID19_GSE168212 [ N = 1 ] , Cell Count = [ 4935 ]</option>
									<option value="COVID19_GSE166992">COVID19_GSE166992 [ N = 9 ] , Cell Count = [ 20348 ]</option>
									<option value="COVID19_GSE154567_severe">COVID19_GSE154567_severe [ N = 3 ] , Cell Count = [ 15190 ]</option>
									<option value="COVID19_GSE154567_moderate">COVID19_GSE154567_moderate [ N = 3 ] , Cell Count = [ 22809 ]</option>
									<option value="COVID19_GSE154567_recovery">COVID19_GSE154567_recovery [ N = 3 ] , Cell Count = [ 29275 ]</option>
									<option value="COVID19_GSE167118_BALF_moderate">COVID19_GSE167118_BALF_moderate [ N = 2 ] , Cell Count = [ 3016 ]</option>
									<option value="COVID19_GSE167118_BALF_severe">COVID19_GSE167118_BALF_severe [ N = 6 ] , Cell Count = [ 23581 ]</option>
									<option value="COVID19_GSE167118_blood_moderate">COVID19_GSE167118_blood_moderate [ N = 1 ] , Cell Count = [ 7270 ]</option>
									<option value="COVID19_GSE167118_blood_severe">COVID19_GSE167118_blood_severe [ N = 6 ] , Cell Count = [ 54094 ]</option>
									<option value="BP_GSE167118_BALF_moderate">BP_GSE167118_BALF_moderate [ N = 3 ] , Cell Count = [ 5560 ]</option>
									<option value="BP_GSE167118_BALF_no">BP_GSE167118_BALF_no [ N = 1 ] , Cell Count = [ 3349 ]</option>
									<option value="BP_GSE167118_blood_moderate">BP_GSE167118_blood_moderate [ N = 3 ] , Cell Count = [ 14866 ]</option>
									<option value="BP_GSE167118_blood_no">BP_GSE167118_blood_no [ N = 1 ] , Cell Count = [ 6988 ]</option>
									<option value="COVID19_GSE169503">COVID19_GSE169503 [ N = 4 ] , Cell Count = [ 4196 ]</option>
									<option value="SLE_GSE142016">SLE_GSE142016 [ N = 3 ] , Cell Count = [ 18531 ]</option>
									<option value="LIHC_GSE107747">LIHC_GSE107747 [ N = 2 ] , Cell Count = [ 9950 ]</option>
									<option value="Lymphoma_GSE124899">Lymphoma_GSE124899 [ N = 2 ] , Cell Count = [ 15124 ]</option>
									<option value="IPEX_GSE167976">IPEX_GSE167976 [ N = 4 ] , Cell Count = [ 10375 ]</option>
					            </select>
							</nav>

            </div>
            
                        <div class="col-md-2 col-sm-2 col-xs-2">
            				<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name">No. of Split</span>
									<span class="menu__item-label"></span>
								</a>
								<select class="language-list-item imm_select" name="top" id="top">
					              <option value="10" selected="selected">10</option>
					              <option value="20">20</option>
					              <option value="30">30</option>
					              <option value="40">40</option>
					              <option value="50">50</option>
					            </select>
							</nav>

            </div>     
            <div class="col-md-3 col-sm-3 col-xs-3">
            				<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name">Score</span>
									<span class="menu__item-label"></span>
								</a>
								<select class="score" name="method" id="method">
					              <option value="normalize" selected="selected">Normalize [ X - min(X) ]</option>
					              <option value="raw">Raw Data</option>
					              <option value="zero">Zero Start [ Xn - X1 ]</option>
					            </select>
							</nav>
            </div>
            <div class="col-md-1 col-sm-1 col-xs-1">
           		<button class="default-btn btn-bg-two border-radius-50 theme-btn" id="submit_botton" type="submit">SUBMIT</button>
            </div>
            </div>
            </form>
		    </section>
		    <section id="tab5" title="Other Pathways">
		    <form action="Code_test.jsp" target="exmdb_function"  method="post" role="form" id="imm_form" name="imm_form">
			<div class="row">
            <div class="col-md-1 col-sm-1 col-xs-1"></div>
            <div class="col-md-10 col-sm-10 col-xs-10">
							<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name"> Other Pathways</span>
									<span class="menu__item-label">(Enter pathway/keywords)</span>
								</a>
								<input type="text" name="searchname_pathway" id="searchname_pathway" class="form-control" value="BIOCARTA_41BB_PATHWAY,BIOCARTA_ACE2_PATHWAY,BIOCARTA_ACH_PATHWAY,BIOCARTA_ACTINY_PATHWAY,BIOCARTA_AGPCR_PATHWAY,BIOCARTA_AGR_PATHWAY,BIOCARTA_AHSP_PATHWAY,BIOCARTA_AKAP13_PATHWAY,BIOCARTA_AKAP95_PATHWAY" placeholder=" NEUTROPHIL_PATHWAY / ACETAMINOPHEN_PATHWAY ..."  style="font-size:15px;padding:20px 15px">
							</nav>
            </div></div>
            <div class="row" style="padding-bottom: 2rem;">
            <div class="col-md-1 col-sm-1 col-xs-1"></div>
            <div class="col-md-5 col-sm-5 col-xs-5">
            				<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name">DATASET ID</span>
									<span class="menu__item-label">(GSE/COVID)</span>
								</a>
								<select class="language-list-item imm_select" name="gseid" id="gseid">
									<option value="BLCA_GSE145281_aPDL1">BLCA_GSE145281_aPDL1 [ N = 10 ] , Cell Count = [ 14474 ]</option>
									<option value="CLL_GSE111014">CLL_GSE111014 [ N = 4 ] , Cell Count = [ 30106 ]</option>
									<option value="KIRC_GSE145281_aPDL1">KIRC_GSE145281_aPDL1 [ N = 4 ] , Cell Count = [ 44220 ]</option>
									<option value="MCC_GSE117988_aPD1aCTLA4">MCC_GSE117988_aPD1aCTLA4 [ N = 1 ] , Cell Count = [ 10148 ]</option>
									<option value="MCC_GSE118056_aPDL1">MCC_GSE118056_aPDL1 [ N = 1 ] , Cell Count = [ 11024 ]</option>
									<option value="NSCLC_GSE127471">NSCLC_GSE127471 [ N = 1 ] , Cell Count = [ 1108 ]</option>
									<option value="PBMC_30K_10X">PBMC_30K_10X [ N = 1 ] , Cell Count = [ 29079 ]</option>
									<option value="PBMC_8K_10X">PBMC_8K_10X [ N = 1 ] , Cell Count = [ 8488 ]</option>
									<option value="PC_GSE67980">PC_GSE67980 [ N = 13 ] , Cell Count = [ 117 ]</option>
									<option value="MEL_GSE157743">MEL_GSE157743 [ N = 22 ] , Cell Count = [ 102 ]</option>
									<option value="COVID19_GSE168212">COVID19_GSE168212 [ N = 1 ] , Cell Count = [ 4935 ]</option>
									<option value="COVID19_GSE166992">COVID19_GSE166992 [ N = 9 ] , Cell Count = [ 20348 ]</option>
									<option value="COVID19_GSE154567_severe">COVID19_GSE154567_severe [ N = 3 ] , Cell Count = [ 15190 ]</option>
									<option value="COVID19_GSE154567_moderate">COVID19_GSE154567_moderate [ N = 3 ] , Cell Count = [ 22809 ]</option>
									<option value="COVID19_GSE154567_recovery">COVID19_GSE154567_recovery [ N = 3 ] , Cell Count = [ 29275 ]</option>
									<option value="COVID19_GSE167118_BALF_moderate">COVID19_GSE167118_BALF_moderate [ N = 2 ] , Cell Count = [ 3016 ]</option>
									<option value="COVID19_GSE167118_BALF_severe">COVID19_GSE167118_BALF_severe [ N = 6 ] , Cell Count = [ 23581 ]</option>
									<option value="COVID19_GSE167118_blood_moderate">COVID19_GSE167118_blood_moderate [ N = 1 ] , Cell Count = [ 7270 ]</option>
									<option value="COVID19_GSE167118_blood_severe">COVID19_GSE167118_blood_severe [ N = 6 ] , Cell Count = [ 54094 ]</option>
									<option value="BP_GSE167118_BALF_moderate">BP_GSE167118_BALF_moderate [ N = 3 ] , Cell Count = [ 5560 ]</option>
									<option value="BP_GSE167118_BALF_no">BP_GSE167118_BALF_no [ N = 1 ] , Cell Count = [ 3349 ]</option>
									<option value="BP_GSE167118_blood_moderate">BP_GSE167118_blood_moderate [ N = 3 ] , Cell Count = [ 14866 ]</option>
									<option value="BP_GSE167118_blood_no">BP_GSE167118_blood_no [ N = 1 ] , Cell Count = [ 6988 ]</option>
									<option value="COVID19_GSE169503">COVID19_GSE169503 [ N = 4 ] , Cell Count = [ 4196 ]</option>
									<option value="SLE_GSE142016">SLE_GSE142016 [ N = 3 ] , Cell Count = [ 18531 ]</option>
									<option value="LIHC_GSE107747">LIHC_GSE107747 [ N = 2 ] , Cell Count = [ 9950 ]</option>
									<option value="Lymphoma_GSE124899">Lymphoma_GSE124899 [ N = 2 ] , Cell Count = [ 15124 ]</option>
									<option value="IPEX_GSE167976">IPEX_GSE167976 [ N = 4 ] , Cell Count = [ 10375 ]</option>
					            </select>
							</nav>

            </div>
            
                        <div class="col-md-2 col-sm-2 col-xs-2">
            				<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name">No. of Split</span>
									<span class="menu__item-label"></span>
								</a>
								<select class="language-list-item imm_select" name="top" id="top">
					              <option value="10" selected="selected">10</option>
					              <option value="20">20</option>
					              <option value="30">30</option>
					              <option value="40">40</option>
					              <option value="50">50</option>
					            </select>
							</nav>

            </div>     
            <div class="col-md-3 col-sm-3 col-xs-3">
            				<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name">Score</span>
									<span class="menu__item-label"></span>
								</a>
								<select class="score" name="method" id="method">
					              <option value="normalize" selected="selected">Normalize [ X - min(X) ]</option>
					              <option value="raw">Raw Data</option>
					              <option value="zero">Zero Start [ Xn - X1 ]</option>
					            </select>
							</nav>
            </div>
            <div class="col-md-1 col-sm-1 col-xs-1">
           		<button class="default-btn btn-bg-two border-radius-50 theme-btn" id="submit_botton" type="submit">SUBMIT</button>     	
            </div>
            </div>
            </form>
		    </section>
		  </div>
        <br>
        <hr>
						<iframe frameborder=0 width="100%" height="450px"
							name="exmdb_function" id="exmdb_function"
							src="Code_test.jsp"
							scrolling="no" onload="load()"></iframe>
  </div>
</div>
<br>
<footer class="footer-area footer-bg">
  <div class="container">
    <div class="copy-right-area">
      <div class="copy-right-text">
        <p>2022,CopyRight © HMU. College of Bioinformatics Science and Technology, Harbin, China.</p>
      </div>
    </div>
  </div>
</footer>
<script src="assets/js/jquery.min.js"></script>
<script src="assets/js/bootstrap.bundle.min.js"></script>
<script src="assets/js/owl.carousel.min.js"></script>
<script src="assets/js/jquery.magnific-popup.min.js"></script>
<script src="assets/js/jquery.nice-select.min.js"></script>
<script src="assets/js/wow.min.js"></script>
<script src="assets/js/meanmenu.js"></script>
<script src="assets/js/jquery.ajaxchimp.min.js"></script>
<script src="assets/js/form-validator.min.js"></script>
<script src="assets/js/contact-form-script.js"></script>
<script src="assets/js/custom.js"></script>






<script type="text/javascript" src="select/demo2.js" ></script>
<script type="text/javascript" src="select/jquery.min.js" ></script>

  <script type="text/javascript" src="tab_search/jquery-tab.js"></script>
	<script type="text/javascript">
		$(function(){
			// Calling the plugin
			$('.tab-group').tabify();
		})
	</script>

<script type="text/javascript" src="select/selectpage.min.js" ></script>
<script type="text/javascript" src="select/data/function_GO_bp.js" ></script>
<script type="text/javascript" src="select/data/function_GO_cc.js" ></script>
<script type="text/javascript" src="select/data/function_GO_mf.js" ></script>
<script type="text/javascript" src="select/data/function_50_hallmaker.js" ></script>
<script type="text/javascript" src="select/data/function_pathway.js" ></script>
<script type="text/javascript">
	$(function(){
		$('#searchname').selectPage({
			showField : 'name',
			keyField : 'id',
			data : tag_data_bp,
			lang: 'en',
			multiple : true,
			maxSelectLimit : 10,
			noResultClean : true
		});
		
		$('#searchname_cc').selectPage({
			showField : 'name',
			keyField : 'id',
			data : tag_data_cc,
			lang: 'en',
			multiple : true,
			maxSelectLimit : 10,
			noResultClean : true
		});
		
		$('#searchname_mf').selectPage({
			showField : 'name',
			keyField : 'id',
			data : tag_data_mf,
			lang: 'en',
			multiple : true,
			maxSelectLimit : 10,
			noResultClean : true
		});
		$('#searchname_hallmaker').selectPage({
			showField : 'name',
			keyField : 'id',
			data : tag_data_hallmaker,
			lang: 'en',
			multiple : true,
			maxSelectLimit : 10,
			noResultClean : true
		});
		
		$('#searchname_pathway').selectPage({
			showField : 'name',
			keyField : 'id',
			data : tag_data_pathway,
			lang: 'en',
			multiple : true,
			maxSelectLimit : 10,
			noResultClean : true
		});
		
		SyntaxHighlighter.all();
	});
</script>



  <script>
    wow = new WOW(
      {
        animateClass: 'animated',
        offset:100
      }
    );
    wow.init();
  </script>
  




</body>
</html>