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

<!-- Topscroll -->
<link rel="stylesheet" href="TopScroll/css/reset.min.css">
<link rel="stylesheet" href="TopScroll/css/style.css">

<link rel="icon" type="image/png" href="assets/images/favicon.png">
<title>ExMdb : Detailed information</title>
<% 
String GSEID = "GSE124158_Bone";
String sql = "";   

if(request.getParameter("GSEID")!=null&&!request.getParameter("GSEID").equals("")){
	GSEID = request.getParameter("GSEID");
}

sql = "select * from data_info where i='"+GSEID+"'";


String searchname = "";
if(request.getParameter("searchname")!=null&&!request.getParameter("searchname").equals("")){	
	searchname = request.getParameter("searchname");
}



Connection conn = new dbcon().getcon();
Statement st = conn.createStatement();
ResultSet rs = null;		
rs = st.executeQuery(sql);
rs.next();

String pubmedid = rs.getString("pubmedid"); 
String gseid = rs.getString("gseid"); 
String mysql_i = rs.getString("i"); 
String dataid = rs.getString("dataid");	 
String organ = rs.getString("organ");
String title = rs.getString("title");
String samplecount = rs.getString("samplecount");
String cancercount = rs.getString("cancercount");
String normalcount = rs.getString("normalcount");
String cancertype = rs.getString("cancertype");
String localtion = rs.getString("localtion");
String rnatype = rs.getString("rnatype");
String dataprocessing = rs.getString("dataprocessing");
String extractionprotocol = rs.getString("extractionprotocol");
String extractionmolecule = rs.getString("extractionmolecule");
String lable = rs.getString("lable");
String labelprotocol = rs.getString("labelprotocol");
String hybridizationprotocol = rs.getString("hybridizationprotocol");
String scanprotocol = rs.getString("scanprotocol");
String datavalue = rs.getString("datavalue");
String summary = rs.getString("summary");
String overalldesign = rs.getString("overalldesign");
String citation = rs.getString("citation");


rs.close();
st.close();
conn.close();

String function_word = "";

function_word = searchname;


String GSVA_show = "display:none";
String GSVA_link = "a.jsp";
String GSVA_title = "";


String imm_show = "display:none";
String imm_link = "a.jsp";
String imm_title = "";

if (rnatype.equals("lncRNA+mRNA")){
	GSVA_show = " ";
	GSVA_link = "exmdb_gsva_50hallmarker.jsp";
	GSVA_title = "<h1 class='module_box_title'>GSVA:50 Hallmarker</h1>";
	imm_show = " ";
	imm_link = "exmdb_imm_cell.jsp";
	imm_title = "<h1 class='module_box_title'>Immune cell infiltration</h1>";
}

%>

<style>

body{
	overflow-x: hidden;
}

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

.imm_select{
	width:100%;
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
      <h3>Detailed information [High-throughput data]</h3>
      <ul>
        <li><a href="index.html">Home</a></li>
        <li><i class='bx bx-chevrons-right'></i></li>
        <li>Detailed information [High-throughput data]</li>
      </ul>
    </div>
  </div>
  <div class="inner-shape"><img src="assets/images/shape/inner-shape.png" alt="Images"></div>
</div>

<div class="terms-conditions-area" style="padding:0px 0px 30px 0px">
<div class="container1">

<h1 class="module_box_title">Data information</h1>
<div style="margin-left: 12.5%;width:25%;margin-bottom:1%;height:5px;background:linear-gradient(to left,#ff000000,#518cffa8,#ff000000);"></div>

<div class="contact-wrapper module_box wow fadeInUp animated">
	<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Data Id :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><span style="color: #657be8;font-weight: 400;"><%=gseid%></span> ( <%=rnatype%> )</div>
		</div>
	</div>
	<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Disease :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><span style="color:#F9B189;font-weight: 400;text-transform:capitalize"><%=cancertype%></span> ( <%=organ%> )</div>
		</div>
	</div>
	<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Title :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><%=title%></div>
		</div>
	</div>
	<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Percentage of samples :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><span style="color: #657be8;font-weight: 400;">Cancer</span> : <%=cancercount%> / <span style="color:#F9B189;font-weight: 400;">Control</span> : <%=normalcount%></div>
		</div>
	</div>
	
	<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Data Processing :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><%=dataprocessing%></div>
		</div>
	</div>

	<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Data Summary :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><%=summary%></div>
		</div>
	</div>
	
		<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Overall design :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><%=overalldesign%></div>
		</div>
	</div>
	
		<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Citation(s) :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><%=citation%></div>
		</div>
	</div>

	<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Extraction Protocol :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><%=extractionprotocol%></div>
		</div>
	</div>

	<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Extraction Molecule :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><%=extractionmolecule%></div>
		</div>
	</div>	


	<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Lable :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><%=lable%></div>
		</div>
	</div>	
	
	<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Label Protocol :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><%=labelprotocol%></div>
		</div>
	</div>
	
	<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Hybridization Protocol :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><%=hybridizationprotocol%></div>
		</div>
	</div>
	
	<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Scan Protocol :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><%=scanprotocol%></div>
		</div>
	</div>
    </div>
<br>

<h1 class="module_box_title">Basic Analysis</h1>
<div style="margin-left: 12.5%;width:25%;margin-bottom:1%;height:5px;background:linear-gradient(to left,#ff000000,#518cffa8,#ff000000);"></div>
<div class="contact-wrapper module_box wow fadeInRight animated">
						<iframe frameborder=0 width="100%" height="1100px"
							name="exmdb_volcano" id="exmdb_volcano"
							src="exmdb_volcano.jsp?GSEID=<%=GSEID%>&searchname=<%=searchname%>"
							scrolling="no" ></iframe></div>

<br>

<h1 class="module_box_title">Functional enrichment</h1>
<div style="margin-left: 12.5%;width:25%;margin-bottom:1%;height:5px;background:linear-gradient(to left,#ff000000,#518cffa8,#ff000000);"></div>
<div class="contact-wrapper module_box wow fadeInRight animated">
						<iframe frameborder=0 width="100%" height="1400px"
							name="exmdb_function" id="exmdb_function"
							src="exmdb_function.jsp?searchname=<%=searchname%>"
							scrolling="no" ></iframe></div>


<br>
<%=GSVA_title%>
<div style="margin-left: 12.5%;width:25%;margin-bottom:1%;height:5px;background:linear-gradient(to left,#ff000000,#518cffa8,#ff000000);<%=GSVA_show%>"></div>
<div class="contact-wrapper module_box wow fadeInRight animated" style="<%=GSVA_show%>">
						<iframe frameborder=0 width="100%" height="800px"
							name="exmdb_function" id="exmdb_function"
							src="<%=GSVA_link%>?cancername=<%=mysql_i%>"
							scrolling="no" ></iframe>
</div>


<br>
<%=imm_title%>
<div style="margin-left: 12.5%;width:25%;margin-bottom:1%;height:5px;background:linear-gradient(to left,#ff000000,#518cffa8,#ff000000);<%=imm_show%>"></div>
<div class="contact-wrapper module_box wow fadeInRight animated" style="<%=imm_show%>">
						<div class="row">
						<div class="col-md-4 col-sm-4 col-xs-12"><h4 style="font-size: 35px;text-align: end;"><span style="font-weight: 400">DATASET : </span><span style="color: #5470c6;font-weight: 900;"><%=gseid%></span></h4></div>
						<div class="col-md-3 col-sm-3 col-xs-12"><h4 style="font-size: 35px;text-align: end;">Method : </h4></div>
						<div class="col-md-5 col-sm-5 col-xs-12">
					            <select class="language-list-item imm_select" name="imm_function" id="imm_function" onchange="changeurl(this.value)">
					              <option value="exmdb_imm_cell.jsp">Cibersort</option>
					              <option value="exmdb_imm_cell_x.jsp">Xcell</option>
					            </select>	
					    </div>
						</div>
						<iframe frameborder=0 width="100%" height="800px"
							name="exmdb_imm" id="exmdb_imm"
							src="<%=imm_link%>?cancername=<%=gseid%>"
							scrolling="no" ></iframe>
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
function changeurl(val)
{
	document.getElementById("exmdb_imm").src=document.getElementById("imm_function").value+"?cancername="+'<%=gseid%>';
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