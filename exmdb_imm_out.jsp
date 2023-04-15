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
<title>ExMdb : Immune cell infiltration</title>


<link rel="stylesheet" type="text/css" href="select/demo.css"/>
<link rel="stylesheet" type="text/css" href="select/style-adsila.css" />
<link rel="stylesheet" href="select/selectpage_cell.css" type="text/css">

<style>
.search_item{
	font-size:20px;
}

.imm_select{
	width:100%;
}


.theme-btn {
	width:100%;
	height:100%;
	padding:2%;
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
  <div class="side-nav-responsive">
    <div class="container-max">
      <div class="dot-menu">
        <div class="circle-inner">
          <div class="in-circle circle-one"></div>
          <div class="in-circle circle-two"></div>
          <div class="in-circle circle-three"></div>
        </div>
      </div>
     
    </div>
  </div>
</div>

<div class="inner-banner">
  <div class="container">
    <div class="inner-title text-center">
      <h3>Immune cell infiltration</h3>
      <ul>
        <li><a href="index.html">Home</a></li>
        <li><i class='bx bx-chevrons-right'></i></li>
        <li>Immune cell infiltration</li>
      </ul>
    </div>
  </div>
  <div class="inner-shape"><img src="assets/images/shape/inner-shape.png" alt="Images"></div>
</div>
<div class="terms-conditions-area" style="padding:40px 0px">
  <div class="container wow fadeInUp animated" style="max-width: 1440px;">
		<div class="contact-header">
			<h1>Immune cell infiltration analysis</h1>
		</div>		
		<div class="contact-wrapper">
	    <form action="exmdb_imm_cell.jsp" target="exmdb_imm"  method="post" role="form" id="imm_form" name="imm_form">
            <div class="row">
            <p class="sub_title_function"><img src="img/iconfont/imm.png" style="margin-top: -7px;margin-right: 7px;">ExMdb-immunity</p>
	    	<p class="sub_text_function">Plot cell infiltration stacks based on sample gene expression. [ Including Cibersort and XCell ].</p>
            <div class="col-md-5 col-sm-5 col-xs-12">
            
            				<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name">Dataset :</span>
									<span class="menu__item-label">( EXORBASE / NCBI-GEO )</span>
								</a>
								<select class="language-list-item imm_select" name="cancername" id="cancername">
					              <option value="HCC">EXORBASE-HCC [ Hepatocellular carcinoma ] ( N = 112 )</option>
					              <option value="BRCA">EXORBASE-BRCA [ Breast cancer ] ( N = 140 )</option>
					              <option value="CRC">EXORBASE-CRC [ Colorectal cancer ] ( N = 35 )</option>
					              <option value="ESCC">EXORBASE-ESCC [ Esophageal squamous cell carcinoma ]( N = 6 )</option>
					              <option value="GBM">EXORBASE-GBM [ Glioblastoma multiforme ] ( N = 13 )</option>
					              <option value="GC">EXORBASE-GC [ Gastric cancer ] ( N = 9 )</option>
								  <option value="KIRC">EXORBASE-KIRC [ Kidney cancer ] ( N = 15 )</option>
								  <option value="MEL">EXORBASE-MEL [ Melanoma ] ( N = 21 )</option>
								  <option value="ML">EXORBASE-ML [ Malignant lymphoma ] ( N = 28 )</option>
								  <option value="OV">EXORBASE-OV [ Ovarian cancer ] ( N = 30 )</option>
								  <option value="PAAD">EXORBASE-PAAD [ Pancreatic adenocarcinoma ] ( N = 164 )</option>
								  <option value="SCLC">EXORBASE-SCLC [ Small cell lung cancer ] ( N = 36 )</option>
								  <option value="GSE125442">GSE125442 [ KIRC ] ( cancer = 10,normal = 10 )</option>
								  <option value="GSE133684">GSE133684 [ Pancreatic ] ( cancer = 284,normal = 117 )</option>
								  <option value="GSE144521">GSE144521 [ Cholangiocarcinoma ] ( cancer = 35,normal = 6 )</option>
								  
								  <option value="GSE93070">GSE93070 [ Breast cancer ] ( cancer = 2,normal = 2 )</option>
								  <option value="GSE106804">GSE106804 [ Glioma ] ( cancer = 8,normal = 6 )</option>
								  <option value="GSE104926">GSE104926 [ Esophageal ] ( cancer = 6,normal = 6 )</option>
					            </select>							
           					 </nav>
            

            </div>
            <div class="col-md-5 col-sm-5 col-xs-12">
                        	<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name">Method :</span>
									<span class="menu__item-label">( Cibersort / Xcell )</span>
								</a>
					            <select class="language-list-item imm_select" name="imm_function" id="imm_function" onchange="changeurl(this.value)">
					              <option value="exmdb_imm_cell.jsp">Cibersort</option>
					              <option value="exmdb_imm_cell_x.jsp">Xcell</option>
					            </select>						
           					 </nav>
            </div>
            <div class="col-md-2 col-sm-2 col-xs-12">
            	<button class="default-btn btn-bg-two border-radius-50 theme-btn" type="submit" style="padding:5px 30px;height:50%;width:100%;font-size:20px;top: 50%;border-style: none;" onclick="show_qi()">Submit</button>
            </div>
            </div>
        </form>
        <br>
		
						<iframe frameborder=0 width="100%" height="800px"
							name="exmdb_imm" id="exmdb_imm"
							src="exmdb_imm_cell.jsp"
							scrolling="no" onload="load()"></iframe></div>
  </div>
</div>
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
function changeurl(val)
{
	document.getElementById("exmdb_imm").src=val;
	document.getElementById("imm_form").action=val;
}
</script>

<script>
function load(){
	$("#graph_qi").fadeOut(1);
}

function show_qi(){
	$("#graph_qi").fadeIn(1);
}
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