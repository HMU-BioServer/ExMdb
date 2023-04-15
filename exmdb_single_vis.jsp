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

<title>ExMdb : Visualization Tools</title>



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

.score{
	width:100%;
}

.score > ul{
	height:200px;
	overflow-y: auto !important;
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
      <h3>Single-Cell Visualization Tools</h3>
      <ul>
        <li><a href="index.html">Home</a></li>
        <li><i class='bx bx-chevrons-right'></i></li>
        <li>Single-Cell Visualization Tools</li>
      </ul>
    </div>
  </div>
  <div class="inner-shape"><img src="assets/images/shape/inner-shape.png" alt="Images"></div>
</div>
  <div class="container wow fadeInUp animated" style="max-width: 1440px;padding-top: 2.5rem;">

		<div class="contact-header">
			<h1>Single-Cell Visualization Tools</h1>
		</div>	
		<div class="contact-wrapper">
			<p class="sub_title_function" style="margin-top:15px"><img src="img/iconfont/vis.png" style="margin-top: -7px;margin-right: 7px;">ExMdb-Single Cell Visualization Tools</p>
	    	<p class="sub_text_function">The single-cell search interface obtains corresponding detailed information by searching for keywords and dataset ID.</p>
				 <nav class="menu menu--adsila">
					<a class="menu__item" href="#">
						<span class="menu__item-name">Seurat</span>
						<span class="menu__item-label">(Integrated Analysis Of ScRNA-seq Datasets)</span>
					</a>
					<div class="row">
						<div class="col-lg-6 col-md-6">
						<label class="label_qi">DataSet ID : </label>
								 <select class="language-list-item score" name="gseid" id="gseid" onchange="change()">
									<option value="BLCA_GSE145281_aPDL1_count">BLCA_GSE145281_aPDL1 [ N = 10 ] , Cell Count = [ 14474 ]</option>
									<option value="CLL_GSE111014_count">CLL_GSE111014 [ N = 4 ] , Cell Count = [ 30106 ]</option>
									<option value="KIRC_GSE145281_aPDL1_count">KIRC_GSE145281_aPDL1 [ N = 4 ] , Cell Count = [ 44220 ]</option>
									<option value="MCC_GSE117988_aPD1aCTLA4_count">MCC_GSE117988_aPD1aCTLA4 [ N = 1 ] , Cell Count = [ 10148 ]</option>
									<option value="MCC_GSE118056_aPDL1_count">MCC_GSE118056_aPDL1 [ N = 1 ] , Cell Count = [ 11024 ]</option>
									<option value="NSCLC_GSE127471_count">NSCLC_GSE127471 [ N = 1 ] , Cell Count = [ 1108 ]</option>
									<option value="PBMC_30K_10X_count">PBMC_30K_10X [ N = 1 ] , Cell Count = [ 29079 ]</option>
									<option value="PBMC_8K_10X_count">PBMC_8K_10X [ N = 1 ] , Cell Count = [ 8488 ]</option>
									<option value="PC_GSE67980_count">PC_GSE67980 [ N = 13 ] , Cell Count = [ 117 ]</option>
									<option value="MEL_GSE157743_count">MEL_GSE157743 [ N = 22 ] , Cell Count = [ 102 ]</option>
									<option value="COVID19_GSE168212_count">COVID19_GSE168212 [ N = 1 ] , Cell Count = [ 4935 ]</option>
									<option value="COVID19_GSE166992_count">COVID19_GSE166992 [ N = 9 ] , Cell Count = [ 20348 ]</option>
									<option value="COVID19_GSE154567_severe_count">COVID19_GSE154567_severe [ N = 3 ] , Cell Count = [ 15190 ]</option>
									<option value="COVID19_GSE154567_moderate_count">COVID19_GSE154567_moderate [ N = 3 ] , Cell Count = [ 22809 ]</option>
									<option value="COVID19_GSE154567_recovery_count">COVID19_GSE154567_recovery [ N = 3 ] , Cell Count = [ 29275 ]</option>
									<option value="COVID19_GSE167118_BALF_moderate_count">COVID19_GSE167118_BALF_moderate [ N = 2 ] , Cell Count = [ 3016 ]</option>
									<option value="COVID19_GSE167118_BALF_severe_count">COVID19_GSE167118_BALF_severe [ N = 6 ] , Cell Count = [ 23581 ]</option>
									<option value="COVID19_GSE167118_blood_moderate_count">COVID19_GSE167118_blood_moderate [ N = 1 ] , Cell Count = [ 7270 ]</option>
									<option value="COVID19_GSE167118_blood_severe_count">COVID19_GSE167118_blood_severe [ N = 6 ] , Cell Count = [ 54094 ]</option>
									<option value="BP_GSE167118_BALF_moderate_count">BP_GSE167118_BALF_moderate [ N = 3 ] , Cell Count = [ 5560 ]</option>
									<option value="BP_GSE167118_BALF_no_count">BP_GSE167118_BALF_no [ N = 1 ] , Cell Count = [ 3349 ]</option>
									<option value="BP_GSE167118_blood_moderate_count">BP_GSE167118_blood_moderate [ N = 3 ] , Cell Count = [ 14866 ]</option>
									<option value="BP_GSE167118_blood_no_count">BP_GSE167118_blood_no [ N = 1 ] , Cell Count = [ 6988 ]</option>
									<option value="COVID19_GSE169503_count">COVID19_GSE169503 [ N = 4 ] , Cell Count = [ 4196 ]</option>
									<option value="SLE_GSE142016_count">SLE_GSE142016 [ N = 3 ] , Cell Count = [ 18531 ]</option>
									<option value="LIHC_GSE107747_count">LIHC_GSE107747 [ N = 2 ] , Cell Count = [ 9950 ]</option>
									<option value="Lymphoma_GSE124899_count">Lymphoma_GSE124899 [ N = 2 ] , Cell Count = [ 15124 ]</option>
									<option value="IPEX_GSE167976_count">IPEX_GSE167976 [ N = 4 ] , Cell Count = [ 10375 ]</option>
					            </select>
						</div>
						<div class="col-lg-6 col-md-6">
						
						<label class="label_qi">Coordinates : </label>
							<select class="language-list-item imm_select" name="location" id="location" onchange="change()">
									<option value="umap" selected="selected">umap</option>
									<option value="tsne">tsne</option>
							</select>
						

						</div>
					</div><hr>
					<div class="row">
						<div class="col-lg-4 col-md-6">
						<label class="label_qi">Resolution : </label>
							<select class="language-list-item imm_select score" name="resolution" id="resolution" onchange="change()">
									<option value="1">0.1</option>
									<option value="2">0.2</option>
									<option value="3">0.3</option>
									<option value="4">0.4</option>
									<option value="5" selected="selected">0.5</option>
									<option value="6">0.6</option>
									<option value="7">0.7</option>
									<option value="8">0.8</option>
									<option value="9">0.9</option>
							</select>
						</div>
						<div class="col-lg-4 col-md-6">
						<label class="label_qi">Group by : </label>
							<select class="language-list-item imm_select" name="character" id="character" onchange="change()">
									<option value="celltype">Cell Type</option>
									<option value="source">Source</option>
									<option value="stage">Stage</option>
									<option value="sample" selected="selected">Sample</option>
									<option value="pseudotime">Pseudotime-Time</option>
									<option value="state">Pseudotime-State</option>
							</select>
						</div>
						<div class="col-lg-4 col-md-6">
						
						<label class="label_qi">Gene Symbol: </label>
							<input type="text" name="flora1" id="flora1" class="form-control" placeholder=" XIST / MALAT1 / ..." onchange="change()">
						

						</div>
					</div>
				</nav>
		<hr>
		
							<iframe frameborder=0 width="100%" height="450px" 
							name="exmdb_cell" id="exmdb_cell"
							src="exmdb_Cell_Map_res.jsp"
							scrolling="no" onload="load()"></iframe>
		
		<hr>
		
		
		<nav class="menu menu--adsila">
					<a class="menu__item" href="#">
						<span class="menu__item-name">Monocle</span>
						<span class="menu__item-label">(Pseudotime Estimation)</span>
					</a>
					<div class="row">
						<div class="col-lg-3 col-md-6"></div>
						<div class="col-lg-6 col-md-6">
						
						<label class="label_qi">DataSet ID : </label>
							   <select class="language-list-item score" name="gseid_time" id="gseid_time" onchange="change_time()">
									<option value="BLCA_GSE145281_aPDL1_count">BLCA_GSE145281_aPDL1 [ N = 10 ] , Cell Count = [ 14474 ]</option>
									<option value="CLL_GSE111014_count">CLL_GSE111014 [ N = 4 ] , Cell Count = [ 30106 ]</option>
									<option value="KIRC_GSE145281_aPDL1_count">KIRC_GSE145281_aPDL1 [ N = 4 ] , Cell Count = [ 44220 ]</option>
									<option value="MCC_GSE117988_aPD1aCTLA4_count">MCC_GSE117988_aPD1aCTLA4 [ N = 1 ] , Cell Count = [ 10148 ]</option>
									<option value="MCC_GSE118056_aPDL1_count">MCC_GSE118056_aPDL1 [ N = 1 ] , Cell Count = [ 11024 ]</option>
									<option value="NSCLC_GSE127471_count">NSCLC_GSE127471 [ N = 1 ] , Cell Count = [ 1108 ]</option>
									<option value="PBMC_30K_10X_count">PBMC_30K_10X [ N = 1 ] , Cell Count = [ 29079 ]</option>
									<option value="PBMC_8K_10X_count">PBMC_8K_10X [ N = 1 ] , Cell Count = [ 8488 ]</option>
									<option value="PC_GSE67980_count">PC_GSE67980 [ N = 13 ] , Cell Count = [ 117 ]</option>
									<option value="MEL_GSE157743_count">MEL_GSE157743 [ N = 22 ] , Cell Count = [ 102 ]</option>
									<option value="COVID19_GSE168212_count">COVID19_GSE168212 [ N = 1 ] , Cell Count = [ 4935 ]</option>
									<option value="COVID19_GSE166992_count">COVID19_GSE166992 [ N = 9 ] , Cell Count = [ 20348 ]</option>
									<option value="COVID19_GSE154567_severe_count">COVID19_GSE154567_severe [ N = 3 ] , Cell Count = [ 15190 ]</option>
									<option value="COVID19_GSE154567_moderate_count">COVID19_GSE154567_moderate [ N = 3 ] , Cell Count = [ 22809 ]</option>
									<option value="COVID19_GSE154567_recovery_count">COVID19_GSE154567_recovery [ N = 3 ] , Cell Count = [ 29275 ]</option>
									<option value="COVID19_GSE167118_BALF_moderate_count">COVID19_GSE167118_BALF_moderate [ N = 2 ] , Cell Count = [ 3016 ]</option>
									<option value="COVID19_GSE167118_BALF_severe_count">COVID19_GSE167118_BALF_severe [ N = 6 ] , Cell Count = [ 23581 ]</option>
									<option value="COVID19_GSE167118_blood_moderate_count">COVID19_GSE167118_blood_moderate [ N = 1 ] , Cell Count = [ 7270 ]</option>
									<option value="COVID19_GSE167118_blood_severe_count">COVID19_GSE167118_blood_severe [ N = 6 ] , Cell Count = [ 54094 ]</option>
									<option value="BP_GSE167118_BALF_moderate_count">BP_GSE167118_BALF_moderate [ N = 3 ] , Cell Count = [ 5560 ]</option>
									<option value="BP_GSE167118_BALF_no_count">BP_GSE167118_BALF_no [ N = 1 ] , Cell Count = [ 3349 ]</option>
									<option value="BP_GSE167118_blood_moderate_count">BP_GSE167118_blood_moderate [ N = 3 ] , Cell Count = [ 14866 ]</option>
									<option value="BP_GSE167118_blood_no_count">BP_GSE167118_blood_no [ N = 1 ] , Cell Count = [ 6988 ]</option>
									<option value="COVID19_GSE169503_count">COVID19_GSE169503 [ N = 4 ] , Cell Count = [ 4196 ]</option>
									<option value="SLE_GSE142016_count">SLE_GSE142016 [ N = 3 ] , Cell Count = [ 18531 ]</option>
									<option value="LIHC_GSE107747_count">LIHC_GSE107747 [ N = 2 ] , Cell Count = [ 9950 ]</option>
									<option value="Lymphoma_GSE124899_count">Lymphoma_GSE124899 [ N = 2 ] , Cell Count = [ 15124 ]</option>
									<option value="IPEX_GSE167976_count">IPEX_GSE167976 [ N = 4 ] , Cell Count = [ 10375 ]</option>
					            </select>
						</div>
					</div><hr>
					
					<div class="row">
					<div class="col-lg-4 col-md-4">
						<label class="label_qi">Resolution : </label>
							<select class="language-list-item imm_select score" name="resolution_2" id="resolution_2" onchange="change_time()">
									<option value="1">0.1</option>
									<option value="2">0.2</option>
									<option value="3">0.3</option>
									<option value="4">0.4</option>
									<option value="5" selected="selected">0.5</option>
									<option value="6">0.6</option>
									<option value="7">0.7</option>
									<option value="8">0.8</option>
									<option value="9">0.9</option>
							</select>
						</div>
						<div class="col-lg-4 col-md-4">
						<label class="label_qi">Group by : </label>
							<select class="language-list-item imm_select" name="character_2" id="character_2" onchange="change_time()">
									<option value="celltype">Cell Type</option>
									<option value="source">Source</option>
									<option value="stage">Stage</option>
									<option value="sample" selected="selected">Sample</option>
									<option value="pseudotime">Pseudotime-Time</option>
									<option value="state">Pseudotime-State</option>
							</select>
						</div>
					
						<div class="col-lg-4 col-md-4">
						<label class="label_qi">Gene Symbol: </label>
							<input type="text" name="time1" id="time1" class="form-control" placeholder=" XIST / MALAT1 / ..." onchange="change_time()">
						</div>
					</div>
					
					
				</nav>
				<hr>
        <br>
						<iframe frameborder=0 width="100%" height="450px"
							name="exmdb_cell_time" id="exmdb_cell_time"
							src="exmdb_Cell_Map_time_res.jsp"
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


<script>
var jq=$.noConflict();
function change()
{
	$("#graph_qi").fadeIn(1);
	document.getElementById("exmdb_cell").src="exmdb_Cell_Map_res.jsp?gseid="+document.getElementById("gseid").value+"&location="+document.getElementById("location").value+"&character="+document.getElementById("character").value+"&resolution="+document.getElementById("resolution").value+"&searchname="+document.getElementById("flora1_text").value;
	

}

function change_time()
{
	$("#graph_qi").fadeIn(1);
	document.getElementById("exmdb_cell_time").src="exmdb_Cell_Map_time_res.jsp?gseid="+document.getElementById("gseid_time").value+"&searchname="+document.getElementById("time1_text").value+"&character="+document.getElementById("character_2").value+"&resolution="+document.getElementById("resolution_2").value;
	
}

function load(){
	$("#graph_qi").fadeOut(1);
}
</script>



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
<script type="text/javascript" src="select/data/highput.js" ></script>
<script type="text/javascript">
	$(function(){
		$('#flora1').selectPage({
			showField : 'name',
			keyField : 'id',
			data : tag_data,
			lang: 'en',
			noResultClean : true
		});
		
		$('#time1').selectPage({
			showField : 'name',
			keyField : 'id',
			data : tag_data,
			lang: 'en',
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