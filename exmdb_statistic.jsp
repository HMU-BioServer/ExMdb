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
<title>ExMdb : Statistic</title>

<!-- Topscroll -->
<link rel="stylesheet" href="TopScroll/css/reset.min.css">
<link rel="stylesheet" href="TopScroll/css/style.css">


<% 
String lncname = "Linc01152";
String pubmedid = "31054511";


String sql = "";   

if(request.getParameter("lncname")!=null&&!request.getParameter("lncname").equals("")){
	lncname = request.getParameter("lncname");
}

if(request.getParameter("pubmedid")!=null&&!request.getParameter("pubmedid").equals("")){
	pubmedid = request.getParameter("pubmedid");
}

sql = "select * from data_biomarker where lncname ='"+lncname+"' and pubmedid ='" + pubmedid + "'";

Connection conn = new dbcon().getcon();
Statement st = conn.createStatement();
ResultSet rs = null;		
rs = st.executeQuery(sql);
rs.next();

String lncname_sql = rs.getString("lncname");
String alias  = rs.getString("alias");
String geneid = rs.getString("geneid");
String lncensg   = rs.getString("lncensg");

String disease = rs.getString("disease");
String cellline = rs.getString("cellline");
String expval = rs.getString("expval");
String regulated = rs.getString("regulated").equals("up-regulated")?"<b>H</b>":"<b>L</b>";

String drug = rs.getString("drug");
String circulating = rs.getString("circulating");
String survival = rs.getString("survival");
String immune = rs.getString("immune");
String metastasis = rs.getString("metastasis");
String recurrence = rs.getString("recurrence");
String cellgrowth = rs.getString("cellgrowth");
String emt = rs.getString("emt");

String apoptosis = rs.getString("apoptosis");
String autophagy = rs.getString("autophagy");
String pubmedid_sql = rs.getString("pubmedid");
String year = rs.getString("year");



rs.close();
st.close();
conn.close();

%>

<style>
.table_contact{
/* 	height:80px; */
	width:95%;
	margin:0rem auto;
	padding:1.5rem;
	font-size: 20px;
	color:#424F60;
    border-bottom: solid 0.5px #9c9c9c;
}

.table_contact:hover{
	box-shadow:0 5px 15px 0 rgb(120 180 214 / 50%);
	transform:translateY(-1px);
	transition:all .5s ease 0s;
	background:aliceblue;
}

.sticky-nav .main-nav {
    -webkit-box-shadow: 0 0 2px rgb(0 0 0);
    box-shadow: 0 0 2px rgb(0 0 0);
}

.module_box{
	width:90%
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


            <li class="nav-item"><a href="exmdb_statistic.jsp" target="_blank" class="nav-link active">Statistic </a></li><li class="nav-item"><a href="exmdb_help.html" target="_blank" class="nav-link">Help </a></li>
            
            <li class="nav-item"><a href="exmdb_download.html" target="_blank" class="nav-link">Download </a></li><li class="nav-item"><a href="exmdb_contact.html" target="_blank" class="nav-link">Contact Us</a></li>
          </ul>

        </div>
      </nav>
    </div>
  </div>

</div>
<div class="search-overlay">
  <div class="d-table">
    <div class="d-table-cell">
      <div class="search-layer"></div>
      <div class="search-layer"></div>
      <div class="search-layer"></div>
      <div class="search-close"><span class="search-close-line"></span><span class="search-close-line"></span></div>
      <div class="search-form">
        <form>
          <input type="text" class="input-search" placeholder="Search here...">
          <button type="submit"><i class='bx bx-search'></i></button>
        </form>
      </div>
    </div>
  </div>
</div>
<div class="inner-banner">
  <div class="container">
    <div class="inner-title text-center">
      <h3>ExMdb - Statistic</h3>
      <ul>
        <li><a href="index.html">Home</a></li>
        <li><i class='bx bx-chevrons-right'></i></li>
        <li>Statistic</li>
      </ul>
    </div>
  </div>
  <div class="inner-shape"><img src="assets/images/shape/inner-shape.png" alt="Images"></div>
</div>



<h1 class="module_box_title" style="margin-left: 5%;">ExMdb Data Composition</h1>
<div style="margin-left: 5%;width:25%;margin-bottom:1%;height:5px;background:linear-gradient(to left,#ff000000,#518cffa8,#ff000000);"></div>
<div class="row contact-wrapper module_box wow fadeInUp animated">
	<div id="container" style="height: 900px;width:60%;"></div>
	<div id="container_bar" style="height: 900px;width:35%;"></div>
</div>
<br><br>
<h1 class="module_box_title" style="margin-left: 5%;">ExMdb human body map</h1>
<div style="margin-left: 5%;width:25%;margin-bottom:1%;height:5px;background:linear-gradient(to left,#ff000000,#518cffa8,#ff000000);"></div>
<div class="contact-wrapper module_box wow fadeInUp animated" style="width:90%">
    	<iframe frameborder=0  name="exmdb_imm" id="exmdb_imm" src="bodymap_sta.jsp" style="width:100%;height:650px" scrolling="no"></iframe>
</div>

<br><br>
<h1 class="module_box_title" style="margin-left: 5%;">Meta data information</h1>
<div style="margin-left: 5%;width:25%;margin-bottom:1%;height:5px;background:linear-gradient(to left,#ff000000,#518cffa8,#ff000000);"></div>
<div class="contact-wrapper module_box wow fadeInUp animated" style="width:90%">
    	<iframe frameborder=0  name="exmdb_imm" id="exmdb_imm" src="exmdb_sta_summary_info_table.jsp" style="width:100%;height:600px" scrolling="no"></iframe>
</div>

<div class="container1">


<br><br>

<h1 class="module_box_title" style="margin-left: 5%;">Sub-cellular localization overview</h1>
<div style="margin-left: 5%;width:25%;margin-bottom:1%;height:5px;background:linear-gradient(to left,#ff000000,#518cffa8,#ff000000);"></div>
<div class="contact-wrapper module_box wow fadeInUp animated">
    	<iframe frameborder=0  name="exmdb_imm" id="exmdb_imm" src="cellmap_all_count.jsp" style="width:90%;margin-left:5%;height:850px" scrolling="no"></iframe>
</div>
<br><br>
</div>

<h1 class="module_box_title" style="margin-left: 5%;">Tools Overview</h1>
<div style="margin-left: 5%;width:25%;margin-bottom:1%;height:5px;background:linear-gradient(to left,#ff000000,#518cffa8,#ff000000);"></div>
<div class="contact-wrapper module_box wow fadeInUp animated" style="width:90%">
    	<iframe frameborder=0  name="exmdb_imm" id="exmdb_imm" src="index_10pic.html" style="width:100%;height:900px" scrolling="no"></iframe>
</div>
<br><br>



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

<script language="javascript">
var $sw = jQuery.noConflict();  
function Change_V(){
	var ind = document.getElementById("GSE_list").value;
	document.getElementById('exmdb_volcano').src = "exmdb_volcano.jsp?GSEID=" + ind;
	document.getElementById('exmdb_function').src = "exmdb_function.jsp?GSEID=" + ind;
}

</script>

<script type="text/javascript" src="echart/echarts.min.js"></script>

<script type="text/javascript">
var dom = document.getElementById("container");
var myChart = echarts.init(dom);
var app = {};

var option;



var data = [
  {
    name: 'High throughput data',
    children: [
      {name:'BRCA',children: [
        {name:'GSE106817 [BRCA,C=115,N=2759]',value: 1},
        {name:'GSE113486 [BRCA,C=40,N=100]',value: 1},
        {name:'GSE118782 [BRCA,C=30,N=10]',value: 1},
        {name:'GSE41922 [BRCA,C=32,N=22]',value: 1},
        {name:'GSE73002 [BRCA,C=1280,N=2686]',value: 1},
        ]},
      {name:'COAD',children: [
        {name:'GSE106817 [CRC,C=115,N=2759]',value: 1},
        {name:'GSE113486 [COAD,C=40,N=100]',value: 1},
        {name:'GSE124158 [COAD,C=30,N=275]',value: 1},
        {name:'GSE149200 [COAD,C=5,N=5]',value: 1},
        {name:'GSE164174 [COAD,C=1423,N=1417]',value: 1},
        {name:'GSE67075 [COAD,C=18,N=30]',value: 1},
        {name:'GSE85589 [COAD,C=5,N=19]',value: 1}
        ]},
      {name:'ESCA',children: [
        {name:'GSE106817 [ESCA,C=88,N=2759]',value: 1},
        {name:'GSE112496 [ESCA,C=5,N=5]',value: 1},
        {name:'GSE112840 [ESCA,C=52,N=52]',value: 1},
        {name:'GSE113486 [ESCA,C=40,N=100]',value: 1},
        {name:'GSE122497 [ESCC,C=566,N=4965]',value: 1},
        {name:'GSE124158 [ESCA,C=30,N=275]',value: 1},
        {name:'GSE164174 [ESCA,C=1423,N=1417]',value: 1},
        {name:'GSE63108 [ESCA,C=28,N=19]',value: 1},
        {name:'GSE104926 [ESCA,C=6,N=6]',value: 1},
        
        ]},
      {name:'STAD',children: [
        {name:'GSE106817 [STAD,C=115,N=2759]',value: 1},
        {name:'GSE113486 [STAD,C=40,N=100]',value: 1},
        {name:'GSE124158 [STAD,C=30,N=275]',value: 1},
        {name:'GSE126399 [STAD,C=12,N=10]',value: 1},
        {name:'GSE130654 [STAD,C=36,N=12]',value: 1},
        {name:'GSE164174 [STAD,C=1423,N=1417]',value: 1},
        {name:'GSE85589 [STAD,C=7,N=19]',value: 1}
      ]},
      {name:'HCC',children: [
        {name:'GSE106817 [HCC,C=115,N=2759]',value: 1},
        {name:'GSE113486 [HCC,C=40,N=100]',value: 1},
        {name:'GSE124158 [HCC,C=30,N=275]',value: 1},
        ]},
      {name:'LUAD',children: [
        {name:'GSE106817 [LUAD,C=115,N=2759]',value: 1},
        {name:'GSE113486 [LUAD,C=40,N=100]',value: 1},
        {name:'GSE124158 [LUAD,C=30,N=275]',value: 1},
        {name:'GSE137140 [LUAD,C=2178,N=1746]',value: 1},
        ]},
      {name:'PAAD',children: [
        {name:'GSE106817 [PAAD,C=115,N=2759]',value: 1},
        {name:'GSE113486 [PAAD,C=40,N=100]',value: 1},
        {name:'GSE124158 [PAAD,C=30,N=275]',value: 1},
        {name:'GSE133684 [PAAD,C=284,N=117]',value: 1},
        {name:'GSE85589 [PAAD,C=88,N=19]',value: 1},
        ]},
      {name:'BSTS',children: [
        {name:'GSE106817 [SARC,C=115,N=2759]',value: 1},
        {name:'GSE110271 [Multiple Myeloma,C=9,N=2]',value: 1},
        {name:'GSE113486 [SARC,C=40,N=100]',value: 1},
        {name:'GSE124158 [BSTS,C=892,N=275]',value: 1},
        ]},
      {name:'CHOL',children: [
        {name:'GSE113486 [CHOL,C=40,N=100]',value: 1},
        {name:'GSE144521 [CHOL,C=35,N=61]',value: 1},
        {name:'GSE85589 [CHOL,C=101,N=19]',value: 1},
        ]},
      {name:'GBM',children: [
        {name:'GSE113486 [GBM,C=40,N=100]',value: 1},
        {name:'GSE122488 [GBM,C=22,N=16]',value: 1},
        {name:'GSE124158 [GBM,C=30,N=275]',value: 1},
        {name:'GSE139031 [GBM,C=423,N=157]',value: 1},
        {name:'GSE106804 [GBM,C=8,N=6]',value: 1},
        ]},
      {name:'PRAD',children: [
        {name:'GSE113486 [PRAD,C=40,N=100]',value: 1},
        {name:'GSE136321 [PRAD,C=24,N=24]',value: 1},
        {name:'GSE138740 [PRAD,C=146,N=89]',value: 1},
        {name:'GSE159177 [PRAD,C=278,N=187]',value: 1},
        ]},
      {name:'OV',children: [
        {name:'GSE106817 [OV,C=320,N=2759]',value: 1},
        {name:'GSE113486 [OV,C=40,N=100]',value: 1},
        {name:'GSE76449 [OV,C=24,N=4]',value: 1},
        ]},
        
      {name:'Othre',children: [
        {name:'GSE113486 [BLCA,C=392,N=100]',value: 1},
        {name:'GSE106817 [SARC,C=115,N=2759]',value: 1},
        {name:'GSE113486 [SARC,C=40,N=100]',value: 1},
        {name:'GSE110271 [Multiple Myeloma,C=9,N=2]',value: 1},
        {name:'GSE125442 [KIRC,C=10,N=10]',value: 1},
        {name:'GSE130512 [Tyriod Cancer,C=16,N=8]',value: 1},
        {name:'GSE37472 [Oral Cancer,C=50,N=26]',value: 1},
        {name:'GSE104926 [Esophagitis,C=6,N=6]',value: 1},
        ]},
    ]
  },
  {
    name: 'Experimental data',
    children: [
      {name: 'Homo sapiens',children: [
        {name:'lung cancer',value: 1},
        {name:'esophageal squamous cell carcinoma',value: 1},
        {name:'bladder cancer',value: 1},
        {name:'colorectal cancer',value: 1},
        {name:'multiple myeloma',value: 1},
        {name:'breast cancer',value: 1},
        {name:'liver cancer',value: 1},
        {name:'pituitary adenoma',value: 1},
        {name:'gastric cancer',value: 1},
        {name:'thyroid cancer',value: 1},
        {name:'pancreatic cancer',value: 1},
        {name:'ovarian cancer',value: 1},
        {name:'glioblastoma',value: 1},
        {name:'oral squamous cell carcinoma',value: 1},
        {name:'cervical cancer',value: 1},
        {name:'prostate cancer',value: 1},
        {name:'colon cancer',value: 1},
        {name:'glioma',value: 1},
        {name:'laryngeal cancer',value: 1},
        {name:'esophageal cancer',value: 1},
        {name:'osteosarcoma',value: 1},
        {name:'cholangiocarcinoma',value: 1},
        {name:'squamous cell carcinoma ',value: 1},
        {name:'Melanoma',value: 1},
        {name:'head and neck squamous cell carcinomas',value: 1},
        {name:'diffuse large b-cell lymphoma',value: 1},
        {name:'intracranial aneurysm ',value: 1},
        {name:'papillary thyroid cancer',value: 1},
        {name:'cholesteatoma',value: 1},
        {name:'endometrial cancer',value: 1},
        {name:'coronary artery aneurysm',value: 1},
        {name:'oral cancer',value: 1},
        {name:'acute myeloid leukemia',value: 1},
        {name:'mantle cell lymphoma',value: 1},
        {name:'chronic myeloid leukemia',value: 1},
        {name:'tongue squamous cell carcinoma',value: 1},
        {name:'urothelial cancer',value: 1},
        {name:'chronic lymphocytic leukemia',value: 1},
        {name:'nasopharyngeal carcinoma',value: 1},
        {name:'renal cancer',value: 1},
        {name:'Advanced colorectal neoplasia',value: 1},
        {name:'Others ... ',value: 1},
        ]
      },
      {name: 'Mus musculus',children: [
        {name:'renal cancer',value: 1},
        {name:'liver cancer',value: 1},
        {name:'skin tumor',value: 1},
        {name:'lung cancer',value: 1},
        {name:'pancreatic cancer',value: 1},
        {name:'Others ... ',value: 1},
        ]
      }
    ]
  },
  {
	    name: 'Single Cell Data',
	    children: [
	      {name: 'Lung disease',children: [
		    {name:'COVID19 GSE168212',value: 1},
		    {name:'COVID19 GSE166992',value: 1},
		    {name:'COVID19 GSE154567 severe',value: 1},
		    {name:'COVID19 GSE154567 moderate',value: 1},
		    {name:'COVID19 GSE154567 recovery',value: 1},
		    {name:'COVID19 GSE167118 BALF moderate',value: 1},
		    {name:'COVID19 GSE167118 BALF severe',value: 1},
		    {name:'COVID19 GSE167118 blood moderate',value: 1},
		    {name:'COVID19 GSE167118 blood severe',value: 1},
		    {name:'BP GSE167118 BALF moderate',value: 1},
		    {name:'BP GSE167118 BALF no',value: 1},
		    {name:'BP GSE167118 blood moderate',value: 1},
		    {name:'BP GSE167118 blood no',value: 1}
	        ]
	      },
	      {name: 'Blood [PBMC]',children: [
			{name:'PBMC 30K 10X',value: 1},
			{name:'PBMC 8K 10X',value: 1},
		        ]
		  },
	      {name: 'Malignant Tumor',children: [
	        {name:'BLCA-GSE145281 aPDL1',value: 1},
	        {name:'CLL-GSE111014',value: 1},
	        {name:'KIRC-GSE145281 aPDL1',value: 1},
	        {name:'MCC-GSE117988 aPD1aCTLA4',value: 1},
	        {name:'MCC-GSE118056 aPDL1',value: 1},
	        {name:'NSCLC-GSE127471',value: 1},
	        {name:'PC-GSE67980',value: 1},
	        {name:'MEL-GSE157743',value: 1},
	        {name:'LIHC-GSE107747',value: 1},
	        {name:'Lymphoma-GSE124899',value: 1},
	        ]
	      },
	      {name: 'Other Diseases',children: [
		    {name:'SLE-GSE142016',value: 1},
		    {name:'IPEX-GSE167976',value: 1},
		    ]
		   }
	    ]
	  },
	  {
		    name: 'Biomarker',
		    children: [
		      {name: 'Drug resistance',children: [{name:'[items = 535]',value: 2},]},
		      {name: 'Circulating',children: [{name:'[items = 878]',value: 2},]},
		      {name: 'Survival',children: [{name:'[items = 4263]',value: 2},]},
		      {name: 'Immune',children: [{name:'[items = 120]',value: 2},]},
		      {name: 'Metastasis',children: [{name:'[items = 3278]',value: 2},]},
		      {name: 'Recurrence',children: [{name:'[items = 292]',value: 2},]},
		      {name: 'Cell Growth',children: [{name:'[items = 2711]',value: 2},]},
		      {name: 'EMT',children: [{name:'[items = 3038]',value: 2},]},
		      {name: 'Apoptosis',children: [{name:'[items = 2829]',value: 2},]},
		      {name: 'Cell Autophagy',children: [{name:'[items = 163]',value: 2},]},
			  
		    ]
		  }
  
];
option = {
		  title: {text: 'Data Statistics of high throughput bulk data, \nsingle-cell data，experimental data and biomarkers data',
				   textStyle: {
					   fontSize: 15,
					   fontFamily: 'Montserrat'
				   }},
  tooltip: {textStyle:{fontFamily: 'Montserrat',fontSize:15}},
  series: {
    type: 'sunburst',
    data: data,
    color: ["#ff5050","#1a53ff","#ffcc00","#ff6600"],
    radius: [0, '90%'],
    itemStyle: {
      borderRadius: 7,
      borderWidth: 2
    },
    levels: [
      {},
      {
        r0: '15%',
        r: '35%',
        itemStyle: {
          borderWidth: 2,
          color: '#ffbf80'
        },
        label: {
          rotate: 'tangential',
          fontFamily: 'Montserrat'
        },
      },
      {
        r0: '35%',
        r: '60%',
        label: {
            fontFamily: 'Montserrat'
          },
      },
     {
        r0: '60%',
        r: '62%',
        label: {
          position: 'outside',
          padding: 3,
          silent: false,
          fontSize: 9,
          fontFamily: 'Montserrat'
        },
        itemStyle: {
          borderWidth: 3
        }
      }

    ]
  }
};

if (option && typeof option === 'object') {
    myChart.setOption(option);
}

</script>

<script type="text/javascript">
var dom = document.getElementById("container_bar");
var myChart = echarts.init(dom);
var app = {};

var option;

option = {
		  title: {text: 'Data Statistics of experimentally validated \nexosome gene-disease associations and biomarkers',
			  
			   textStyle: {
				   fontSize: 15,
				   fontFamily: 'Montserrat'
			   }},
  dataset: {
    source: [
      ['score', 'amount', 'product'],
      [5684,5684,'Others .. '],
      [28,28,'Osteoarthritis'],
      [28,28,'Liver Fibrosis'],
      [36,36,'Alzheimer Disease'],
      [68,68,'multiple myeloma'],
      [89,89,'cholangiocarcinoma'],
      [114,114,'Liver Cancer'],
      [169,169,'Renal Cancer'],
      [202,202,'Papillary thyroid cancer'],
      [236,236,'Squamous Cell Carcinoma'],
      [266,266,'Esophageal Cancer'],
      [359,359,'Bladder Cancer'],
      [383,383,'Cervical Cancer'],
      [384,384,'Osteosarcoma'],
      [391,391,'Ovarian Cancer'],
      [483,483,'Glioblastoma'],
      [543,543,'Malignant Melanoma'],
      [586,586,'Colon Cancer'],
      [659,659,'Prostate Cancer'],
      [1136,1136,'Breast Cancer'],
      [1137,1137,'Gastric Cancer'],
      [1234,1234,'Pancreatic adenocarcinoma '],
      [1375,1375,'Lung Cancer'],
      [1857,1857,'Colorectal Cancer'],
      [2123,2123,'Hepatocellular Carcinoma'],
      

    ]
  },
  tooltip: {textStyle:{fontFamily: 'Montserrat',fontSize:15},},
  grid: { containLabel: true },
  xAxis: { name: 'count',axisLabel: {fontFamily: "Montserrat"},},
  yAxis: { type: 'category',axisLabel: {fontFamily: "Montserrat"},},
  visualMap: {
    orient: 'horizontal',
    left: 'center',
    min: 1,
    max: 1023,
    textStyle:{fontFamily: 'Montserrat',fontSize:12},
    text: [' ', 'Counts'],
    // Map the score column to color
    dimension: 0,
    inRange: {
        color: ['#F5D043','#114182']
      }
  },
  series: [
    {
      type: 'bar',
        label: {
          show: true,
          precision: 2,
          position: 'right',
          valueAnimation: true,
          fontFamily: 'Montserrat'
        },
      encode: {
        // Map the "amount" column to X axis.
        x: 'score',
        // Map the "product" column to Y axis
        y: 'product'
      }
    }
  ]
};

if (option && typeof option === 'object') {
    myChart.setOption(option);
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

<!-- Topscroll -->
<script src="js/1.10.2/jquery.min.js"></script>
<script src="TopScroll/js/index.js"></script>

</body>
</html>