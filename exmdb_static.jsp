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
<title>ExMdb : Data information [Biomarker]</title>

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
            <li class="nav-item"><a href="index.html" class="nav-link active">Home</a></li>
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
                <li class="nav-item"><a href="exmdb_net_out.jsp" target="_blank" class="nav-link">ExMdb-Network </a></li><li class="nav-item"><a href="exmdb_blast_out.jsp" target="_blank" class="nav-link">ExMdb-BLAST </a></li><li class="nav-item"><a href="exmdb_function_time.jsp" target="_blank" class="nav-link">ExMdb-Pseudotime pathway</a></li>
				<li class="nav-item"><a href="#" class="nav-link">ExMdb-FunctionEnrichment <i class='bx bx-caret-down'></i></a><ul class="dropdown-menu"><li class="nav-item"><a href="exmdb_function_mRNA_out.jsp" target="_blank" class="nav-link" style="text-transform: none;">mRNA-FunctionEnrichment</a></li><li class="nav-item"><a href="exmdb_function_out_miRNA.jsp" target="_blank" class="nav-link" style="text-transform: none;">miRNA-FunctionEnrichment</a></li> <li class="nav-item"><a href="exmdb_function_out.jsp" target="_blank" class="nav-link" style="text-transform: none;">lncRNA-FunctionEnrichment</a></li></ul></li>
              </ul>
            </li>


            <li class="nav-item"><a href="exmdb_help.html" target="_blank" class="nav-link">Help </a></li>
            
            <li class="nav-item"><a href="exmdb_contact.html" target="_blank" class="nav-link">Contact Us</a></li>
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
      <h3>Data information [Biomarker]</h3>
      <ul>
        <li><a href="index.html">Home</a></li>
        <li><i class='bx bx-chevrons-right'></i></li>
        <li>Data information [Biomarker]</li>
      </ul>
    </div>
  </div>
  <div class="inner-shape"><img src="assets/images/shape/inner-shape.png" alt="Images"></div>
</div>




<div class="container1">


<h1 class="module_box_title">Data information [ Biomarker ]</h1>
<div style="margin-left: 12.5%;width:25%;margin-bottom:1%;height:5px;background:linear-gradient(to left,#ff000000,#518cffa8,#ff000000);"></div>


<div class="contact-wrapper module_box wow fadeInUp animated">
	<div class="row">
	<div class="col-xl-7 col-lg-7 col-md-7 col-sm-7">
	<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Gene name :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><span style="color: #657be8;font-weight: 400;text-transform:uppercase"><%=lncname_sql%> [ <%=lncensg%> ] [ <%=alias%> ]</span></div>
		</div>
	</div>
	<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Disease :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><span style="color:#ecb71f;font-weight: 400;text-transform:capitalize"><%=disease%> [ <%=cellline%> ]</span></div>
		</div>
	</div>

	<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box" style="font-size:18px">Experimental method :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><span style="color:#4d4dff;font-weight: 400;"><%=expval%></span></div>
		</div>
	</div>
	
	<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Regulated :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><span><%=regulated%></span></div>
		</div>
	</div>
	
	<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Pubmed ID :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><%=pubmedid_sql%> [ <%=year%> ]</div>
		</div>
	</div>
	
	</div>
	
	<div class="col-xl-5 col-lg-5 col-md-5 col-sm-5">
		<div id="container" style="height: 100%;border-left: 1px solid darkgrey;"></div>
	</div>
	</div>
	
	
 </div>

<br><br>


<h1 class="module_box_title">Other links</h1>
<div style="margin-left: 12.5%;width:25%;margin-bottom:1%;height:5px;background:linear-gradient(to left,#ff000000,#518cffa8,#ff000000);"></div>

<div class="contact-wrapper module_box wow fadeInUp animated">
    	<h4 class="home_title" style="color:#424F60">External Annotation for <span style="color:#F9B189;text-transform:uppercase"><%=lncname_sql%></span></h4><hr>
    	<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box"><a href="http://www.noncode.org/gene_trans_search.php?keyword=<%=lncname_sql%>"><img class='table_img' style='width:200px;' src="img/links/noncode.png"></a></div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box">An integrated knowledge database dedicated to ncRNAs, especially lncRNAs.</div>
		</div>
		</div>
		
		<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box" ><a href="https://pubmed.ncbi.nlm.nih.gov/?term=<%=lncname_sql%>"><img class='table_img' style='width:190px;' src="img/links/pubmed.png"></a></div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box">PubMed® comprises more than 32 million citations for biomedical literature from MEDLINE, life science journals, and online books.</div>
		</div>
		</div>
		
		
		<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box"><a href="https://www.ncbi.nlm.nih.gov/genome/?term=<%=lncname_sql%>"><img class='table_img' style='width:190px;' src="img/links/NCBI.png"></a></div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box">National Center for Biotechnology Information.</div>
		</div>
		</div>

		<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box"><a href="https://www.gencodegenes.org/"><img class='table_img' style='width:190px;' src="img/links/Gencode.png"></a> 	</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box">GENCODE are supporting the annotation of non-canonical human ORFs predicted by Ribo-seq data.</div>
		</div>
		</div>

		<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box" ><a href="https://www.genecards.org/lookup/<%=lncname_sql%>"><img class='table_img' style='width:190px;' src="img/links/genecard.png"></a> 	</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box" >GeneCards enables researchers to effectively navigate and inter-relate the wide universe of human genes, diseases etc.</div>
		</div>
		</div>

		<h4 class="home_title" style="color:#424F60;margin-top:4rem">External Annotation for <span style="color:#657be8;text-transform:capitalize"><%=disease%></span></h4><hr>
		<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box"><a href="https://pubmed.ncbi.nlm.nih.gov/?term=<%=disease%>"><img class='table_img' style='width:190px;' src="img/links/pubmed.png"></a></div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box">PubMed® comprises more than 32 million citations for biomedical literature from MEDLINE, life science journals, and online books.</div>
		</div>
		</div>

		<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box"><a href="https://omim.org/search?search=<%=disease%>"><img class='table_img' style='width:190px;' src="img/links/OMIM.png"></a></div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box">A Knowledgebase of Human Genes and Genetic Phenotypes</div>
		</div>
		</div>
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



option = {
  title: {
	   text: ' Hallmarker Radar Chart',
	   textStyle: {
		   fontSize: 15,
		   fontFamily: 'Montserrat'
	   }
   },
   tooltip: {
       show: true,
       trigger: "item",
       textStyle:{fontFamily: 'Montserrat',fontSize:15},
   },
   toolbox: {
       show : true,
       feature : {
           mark : {show: true},
           dataView : {show: true, readOnly: false},
           saveAsImage : {show: true}
       }
   },
  radar: {
    shape: 'circle',
    axisName: {
        color: 'black',
        fontSize: 14,
        fontFamily: 'Montserrat'
      },
    indicator: [
      { name: 'Drug', max: 1,min:-1 },
      { name: 'Circulating', max: 1,min:-1 },
      { name: 'Survival', max:1,min:-1 },
      { name: 'Immune', max: 1,min:-1 },
      { name: 'Metastasis', max: 1,min:-1 },
      { name: 'Recurrence', max: 1,min:-1 },
      { name: 'Cellgrowth', max:1,min:-1 },
      { name: 'Emt', max: 1,min:-1},
      { name: 'Apoptosis', max: 1,min:-1 },
      { name: 'Autophagy', max: 1,min:-1}
    ]
  },
  series: [
    {
      name: 'Hallmarker',
      type: 'radar',
      data: [
        {
          value: [<%=drug%>, <%=circulating%>, <%=survival%>,<%=immune%>,<%=metastasis%>,<%=recurrence%>,<%=cellgrowth%>,<%=emt%>,<%=apoptosis%>,<%=autophagy%>,],
          name: 'Current Gene',
          areaStyle: {
              normal: {
                  color: "#4dadff6b",
              }
          },
          itemStyle: {
              color: "#003d99",
              borderColor:'#4dadff6b',
              borderWidth:10,
          },
          symbolSize: 10,
        },
      ]
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