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
<link rel="stylesheet" href="select/selectpage.css" type="text/css">

<link rel="stylesheet" type="text/css" href="bootstrap_select/bootstrap.min.css">
<link rel="stylesheet" type="text/css" href="bootstrap_select/bootstrap-select.css">

<title>ExMdb : Survival Analyze Tool</title>



<style>
.imm_select{
	width:100%;
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

<% 
String searchname = "HOTAIR";  
if(request.getParameter("searchname")!=null&&!request.getParameter("searchname").equals("")){
	searchname = request.getParameter("searchname");
}

String organ = "lung";  
if(request.getParameter("organ")!=null&&!request.getParameter("organ").equals("")){
	organ = request.getParameter("organ");
}


// String sql = "select * from exmdb_net where node1 in (\"A1BG-AS1\",\"MALAT1\",\"XIST\",\"miR-3605-5p\") or node2 in (\"A1BG-AS1\",\"MALAT1\",\"XIST\",\"miR-3605-5p\")";
// String data = dbhello.exmdb_netmaker(sql);
// System.out.println(data);
%>


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
            <li class="nav-item"><a href="#" class="nav-link">Tools <i class='bx bx-caret-down'></i></a>
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
      <h3>Survival Analyze Tool</h3>
      <ul>
        <li><a href="index.html">Home</a></li>
        <li><i class='bx bx-chevrons-right'></i></li>
        <li>Survival Analyze Tool</li>
      </ul>
    </div>
  </div>
  <div class="inner-shape"><img src="assets/images/shape/inner-shape.png" alt="Images"></div>
</div>
<div class="terms-conditions-area">
<div class="container wow fadeInUp animated" style="max-width: 1440px;padding:30px 0px">
		<div class="contact-header">
			<h1>Survival Analyze Tool</h1>
		</div>	
		<div class="contact-wrapper" style="padding:30px">
<!-- 		<form id="gene" name="gene" method="post" action="http://localhost:8080/Server_survival/CellTracker_Survival_server.jsp" target="exmdb_net"> -->
	<form id="gene" name="gene" method="post" action="http://www.bio-server.cn/Server_survival/CellTracker_Survival_server.jsp" target="exmdb_net">
		
		    <p class="sub_title_function"><img src="img/iconfont/survival.png" style="margin-top: -7px;margin-right: 7px;width: 27px;">ExMdb-Survival</p>
	    	<p class="sub_text_function">A tool to perform COX regression analysis and survival curves for RNA across more than 30 types of malignant cancers.</p>
				 <nav class="menu menu--adsila">
					<a class="menu__item" href="#">
						<span class="menu__item-name">Survival Analyze Tool</span>
						<span class="menu__item-label">(Gene/Data Set/Group color)</span>
					</a>
					
					<div class="row">
						<div class="col-lg-3 col-md-3">
						<label>Gene (Symbol) : </label>
							 <input type="text" name="genename" id="genename" class="form-control"  >
						</div>

						<div class="col-lg-3 col-md-3">
						<label>Data Set : </label>
						<select class="selectpicker show-menu-arrow" style="width:100%" data-live-search="true" name="dis" id="dis" data-size="5">
						<optgroup label="Adrenocortical Carcinoma">
						<option value="ACC">ACC-TCGA (N=79)</option> </optgroup>
						<optgroup label="Bladder Urothelial Carcinoma">
						<option value="BLCA">BLCA-TCGA (N=404)</option> 
						<option value="GSE13507BLCA">BLCA-GSE13507 (N=165)</option>
						<option value="GSE31684BLCA">BLCA-GSE31684 (N=77)</option>
						</optgroup>
							<optgroup label="Breast Invasive Carcinoma">
									<option value="BRCA">BRCA-TCGA (N=1089)</option> 
									<option value="GSE2603BRCA">BRCA-GSE2603 (N=82)</option> 
									<option value="GSE5327BRCA">BRCA-GSE5327 (N=58)</option> 
									<option value="GSE6130BRCA">BRCA-GSE6130 (N=107)</option> 
									<option value="GSE10886BRCA">BRCA-GSE10886 (N=101)</option> 
									<option value="GSE10893BRCA">BRCA-GSE10893 (N=237)</option> 
									<option value="GSE18229BRCA">BRCA-GSE18229 (N=253)</option> 
									<option value="GSE20685BRCA">BRCA-GSE20685 (N=327)</option> 
									<option value="GSE20713BRCA">BRCA-GSE20713 (N=88)</option>
									<option value="GSE21653BRCA">BRCA-GSE21653 (N=252)</option>
									<option value="GSE31448BRCA">BRCA-GSE31448 (N=246)</option>
									<option value="GSE37181BRCA">BRCA-GSE37181 (N=123)</option>
									<option value="GSE42568BRCA">BRCA-GSE42568 (N=104)</option>
									<option value="GSE48390BRCA">BRCA-GSE48390 (N=81)</option>
									<option value="GSE58812BRCA">BRCA-GSE58812 (N=107)</option>
									<option value="GSE88770BRCA">BRCA-GSE88770 (N=117)</option>
							</optgroup>
							<optgroup label="Cervical Squamous Cell Carcinoma">
									<option value="CESC">CESC-TCGA (N=304)</option> </optgroup>
									<optgroup label="Cholangiocarcinoma">
									<option value="CHOL">CHOL-TCGA (N=36)</option> </optgroup>
									<optgroup label="Colon Adenocarcinoma">
									<option value="COAD">COAD-TCGA (N=454)</option> 
									<option value="GSE38882COAD">COAD-GSE38882 (N=122)</option> 
									
									</optgroup>
									<optgroup label="Bladder Urothelial Carcinoma">
									<option value="DLBC">DLBC-TCGA (N=47)</option></optgroup>
							<optgroup label="Esophageal Carcinoma">
									<option value="ESCA">ESCA-TCGA (N=163)</option>
									<option value="GSE53622ESCA">ESCA-GSE53622 (N=60)</option> 
									<option value="GSE53624ESCA">ESCA-GSE53624 (N=119)</option> 
							</optgroup>
							<optgroup label="Glioblastoma Multiforme">
									<option value="GBM">GBM-TCGA (N=159)</option> </optgroup>
									<optgroup label="Head and Neck Squamous Cell Carcinoma">
									<option value="HNSC">HNSC-TCGA (N=500)</option> </optgroup>
									<optgroup label="Kidney Chromophobe">
									<option value="KICH">KICH-TCGA (N=65)</option> </optgroup>
									<optgroup label="Kidney Renal Clear Cell Carcinoma">
									<option value="KIRC">KIRC-TCGA (N=529)</option> </optgroup>
									<optgroup label="Kidney Renal Papillary Cell Carcinoma">
									<option value="KIRP">KIRP-TCGA (N=289)</option></optgroup>
									<optgroup label="Acute Myeloid Leukemia">
									<option value="LAML">LAML-TCGA (N=151)</option> </optgroup>
									<optgroup label="Brain Lower Grade Glioma">
									<option value="LGG">LGG-TCGA (N=506)</option> </optgroup>
									<optgroup label="Liver Hepatocellular Carcinoma">
									<option value="LIHC">LIHC-TCGA (N=370)</option> </optgroup>
							<optgroup label="Lung Adenocarcinoma">
									<option value="LUAD">LUAD-TCGA (N=515)</option> 
									<option value="GSE3141Lung">LUAD-GSE3141 (N=111)</option>
									<option value="GSE11969Lung">LUAD-GSE11969 (N=149)</option> 
									<option value="GSE31210Lung">LUAD-GSE31210 (N=226)</option> 
									<option value="GSE26939Lung">LUAD-GSE26939 (N=115)</option> 
									<option value="GSE30219Lung">LUAD-GSE30219 (N=293)</option> 
							</optgroup>
									
							<optgroup label="Lung Squamous Cell Carcinoma">		
									<option value="LUSC">LUSC-TCGA(N=497)</option> 
									<option value="GSE8894Lung">LUSC-GSE8894 (N=138)</option> 
									<option value="GSE19188Lung">LUSC-GSE19188 (N=82)</option> 
									<option value="GSE37745Lung">LUSC-GSE37745 (N=196)</option> 
									<option value="GSE50081Lung">LUSC-GSE50081 (N=181)</option> 
									<option value="GSE42127Lung">LUSC-GSE42127 (N=176)</option> 							
									
							</optgroup>
							<optgroup label="Mesothelioma">
									<option value="MESO">MESO-TCGA (N=85)</option> </optgroup>
									<optgroup label="Ovarian Serous Cystadenocarcinoma">
									<option value="OV">OV-TCGA (N=375)</option> 
									<option value="GSE14764OV">OV-GSE14764 (N=80)</option> 
									<option value="GSE17260OV">OV-GSE17260 (N=110)</option> 
									<option value="GSE26712OV">OV-GSE26712 (N=56)</option> 
									<option value="GSE30161OV">OV-GSE30161 (N=58)</option> 
									<option value="GSE31245OV">OV-GSE31245 (N=57)</option> 
									<option value="GSE32602OV">OV-GSE32602 (N=260)</option> 
									<option value="GSE49997OV">OV-GSE49997 (N=194)</option> 
									<option value="GSE63885OV">OV-GSE63885 (N=70)</option> 
							</optgroup>
									<optgroup label="Pancreatic Adenocarcinoma">
									<option value="PAAD">PAAD-TCGA (N=177)</option> </optgroup>
									<optgroup label="Pheochromocytoma and Paraganglioma">
									<option value="PCPG">PCPG-TCGA (N=179)</option> </optgroup>
									<optgroup label="Prostate Adenocarcinoma">
									<option value="PRAD">PRAD-TCGA (N=475)</option> 
									<option value="GSE16560PRAD">PRAD-GSE16560 (N=281)</option> 
									</optgroup>
									<optgroup label="Rectum Adenocarcinoma">
									<option value="READ">READ-TCGA (N=166)</option> </optgroup>
									<optgroup label="Sarcoma">
									<option value="SARC">SARC-TCGA (N=257)</option> </optgroup>
									<optgroup label="Stomach Adenocarcinoma">
									<option value="SKCM">SKCM-TCGA (N=456)</option> </optgroup>
							<optgroup label="Stomach Adenocarcinoma">
									<option value="STAD">STAD-TCGA (N=380)</option> 
									<option value="GSE15459STAD">STAD-GSE15459 (N=192)</option> 
									<option value="GSE62254STAD">STAD-GSE62254 (N=300)</option> 
							</optgroup>
							<optgroup label="Testicular Germ Cell Tumors">
									<option value="TGCT">TGCT-TCGA (N=133)</option> 
							</optgroup>
							<optgroup label="Thyroid Carcinoma">
									<option value="THCA">THCA-TCGA (N=502)</option> 
							</optgroup>
							<optgroup label="Thymoma">
									<option value="THYM">THYM-TCGA (N=119)</option> </optgroup>
							<optgroup label="Uterine Carcinosarcoma">
									<option value="UCEC">UCEC-TCGA (N=542)</option> </optgroup>
									<optgroup label="Uterine Corpus Endometrial Carcinoma">
									<option value="UCS">UCS-TCGA (N=56)</option> </optgroup>
									<optgroup label="Uveal Melanoma">
									<option value="UVM">UVM-TCGA (N=80)</option> </optgroup>					 
									</select>
						</div>

						<div class="col-lg-2 col-md-6">
						<label>Egde width : </label>
							<select class="selectpicker show-menu-arrow" id="linel" name="linel">
							<optgroup label="WIDE">
							<option value="6">[Wide] 6</option> 
							<option value="5">[Wide] 5</option> 
							</optgroup>
							<optgroup label="MEDIUM">
							<option value="4">[Medium] 4</option> 
							<option value="3" selected="selected">[Medium] 3</option>
							</optgroup>
							<optgroup label="NARROW">
							<option value="2">[Narrow] 2</option> 
							<option value="1">[Narrow] 1</option>
							</optgroup>
							</select>	
						</div>
						
						<div class="col-md-1 col-sm-1 col-xs-1">
						<label style="font-size: 12px;">High Risk: </label>
							<div class="imm_select">
								<input id="highc" name="highc" type="color" value="#e60000" style="height:35px;width:100%">  
							</div>
						</div>
						
						<div class="col-md-1 col-sm-1 col-xs-1">
						<label style="font-size: 12px;">Low Risk: </label>
							<div class="imm_select">
								<input name="lowc" id="lowc" type="color" value="#002db3" style="height:35px;width:100%">
							</div>
						</div>
						
						<div class="col-lg-2 col-md-2">
							<div class="progress-button elastic">
									<button class="default-btn btn-bg-two border-radius-50" type="submit" style="width: 100%;margin-top: 7.5%;" onclick="show_qi()"><span>Submit</span></button>
							</div>
						</div>
					</div>
				</nav>
		</form>
		<hr>
							<iframe frameborder=0 width="100%" height="600px" 
							name="exmdb_net" id="exmdb_net"
							src="http://www.bio-server.cn/Server_survival/CellTracker_Survival_server.jsp?genename=EGR1&highc=red&lowc=blue&linel=1"
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
<script src="assets/js/noc_jquery.nice-select.min.js"></script>
<script src="assets/js/wow.min.js"></script>
<script src="assets/js/meanmenu.js"></script>
<script src="assets/js/jquery.ajaxchimp.min.js"></script>
<script src="assets/js/form-validator.min.js"></script>
<script src="assets/js/contact-form-script.js"></script>
<script src="assets/js/custom.js"></script>


<script src="http://cdn.bootcss.com/jquery/1.11.0/jquery.min.js" type="text/javascript"></script>
<script>window.jQuery || document.write('<script src="js/jquery-1.11.0.min.js"><\/script>')</script>
<script src="http://cdn.bootcss.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>
<script src="bootstrap_select/bootstrap-select.js"></script>


<script type="text/javascript" src="select/demo2.js" ></script>
<script type="text/javascript" src="select/jquery.min.js" ></script>
<script type="text/javascript" src="select/selectpage.min.js" ></script>
<script type="text/javascript" src="select/data/lncmRNA_survival.js" ></script>


<script type="text/javascript">
	$(function(){
		$('#genename').selectPage({
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

<script>
function load(){
	$("#graph_qi").fadeOut(1);
}

function show_qi(){
	$("#graph_qi").fadeIn(1);
}
</script>



</body>
</html>