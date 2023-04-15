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
<link rel="stylesheet" href="font-awesome-4.7.0/css/font-awesome.css">
<link rel="icon" type="image/png" href="assets/images/favicon.png">
<title>Techex - Technology Services HTML Template</title>


<style>
.title_box{
	font-weight: 700;
}

body::-webkit-scrollbar-track
{
/* 	-webkit-box-shadow: inset 0 0 6px rgba(0,0,0,0.1); */
	background-color: white;
	border-radius: 10px;
}

body::-webkit-scrollbar
{
	width: 5px;
	background-color: white;
}

body::-webkit-scrollbar-thumb
{
	border-radius: 10px;
	width: 5px;
	height: 5px;
	background-color: white;
}

.in_span{
	display: contents !important;
	text-transform: capitalize;
}

.module_box{
	width:100%;
}
</style>

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

.Deatil_link{
	margin-left: 1rem;
	font-size: 15px;
    background: #19359a;
    border: white solid 2px;
    color: white;
    padding: 0.5rem;
    border-radius: 10px;
    cursor: pointer;
}

.Deatil_link:hover{
	transition:all .5s ease 0s;
	background:white;
	border: #19359a solid 2px;
}

.box_link{
	display: block;
	color: inherit;
}
</style>

<% 




String searchname = "BLCA_GSE145281_aPDL1";
if(request.getParameter("searchname")!=null&&!request.getParameter("searchname").equals("")){	
	searchname = request.getParameter("searchname");
}

String sql = "";   
sql = "select * from sc_data_info where datasetname ='"+searchname+"'";

Connection conn = new dbcon().getcon();
Statement st = conn.createStatement();
ResultSet rs = null;		
rs = st.executeQuery(sql);
rs.next();

String datasetname = rs.getString("datasetname"); 
String diseasename1 = rs.getString("diseasename1");
String diseasename2 = rs.getString("diseasename2");
String tcga1 = rs.getString("tcga1");
String tcga2 = rs.getString("tcga2");
String species = rs.getString("species");
String organ_name = rs.getString("organ");
String accession = rs.getString("accession");
String treatment = rs.getString("treatment");
String patients = rs.getString("patients");
String cellcount = rs.getString("cellcount");
String platform = rs.getString("platform");
String primeta = rs.getString("primeta");
String publication = rs.getString("publication");
String pubmedid = rs.getString("pubmedid");
String proteincodingcounts = rs.getString("proteincodingcounts");
String lncrnacounts = rs.getString("lncrnacounts");
String datasetnamelong = rs.getString("datasetnamelong");


rs.close();
st.close();
conn.close();

// System.out.println("Up: "+ up_count);
// System.out.println("Down: "+ down_count);

%>

</head>
<body>
<div class="preloader">
  <div class="d-table">
    <div class="d-table-cell">
      <div class="spinner"></div>
    </div>
  </div>
</div>




<section class="services-style-area pt-50 pb-50">
<div class="container">

<h1 style="margin-left: -5%;" class="module_box_title"><span>Data information</span><a href="exmdb_sc_detail.jsp?gseid=<%=datasetnamelong%>" class="Deatil_link" target="_blank">GET Detail  <i class="fa fa-send"></i></a></h1>
<div style="width:25%;margin-bottom:1%;height:5px;background:linear-gradient(to left,#ff000000,#518cffa8,#ff000000);"></div>


<div class="contact-wrapper module_box wow fadeInUp animated">
	<div class="table_contact">
	<a class="box_link" href="exmdb_sc_detail.jsp?gseid=<%=datasetnamelong%>" target="_blank">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Data Id :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><span style="color: #657be8;font-weight: 400;"><%=datasetname%></span> ( <%=accession%> )</div>
		</div></a>
	</div>
	<div class="table_contact">
	<a class="box_link" href="exmdb_sc_detail.jsp?gseid=<%=datasetnamelong%>" target="_blank">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Disease :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><span style="color:#F9B189;font-weight: 400;text-transform:capitalize"><%=diseasename1%></span> ( <%=diseasename2%> )</div>
		</div></a>
	</div>
	<div class="table_contact">
	<a class="box_link" href="exmdb_sc_detail.jsp?gseid=<%=datasetnamelong%>" target="_blank">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">TCGA :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><span style="color:#F9B189;font-weight: 400;text-transform:capitalize"><%=tcga1%></span> ( <%=tcga2%> )</div>
		</div></a>
	</div>
	<div class="table_contact">
	<a class="box_link" href="exmdb_sc_detail.jsp?gseid=<%=datasetnamelong%>" target="_blank">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Species  :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><span style="color: #657be8;font-weight: 400;"><%=species%></span> ( N = <%=patients%> )</div>
		</div></a>
	</div>
	<div class="table_contact">
	<a class="box_link" href="exmdb_sc_detail.jsp?gseid=<%=datasetnamelong%>" target="_blank">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Gene composition :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><span style="color: #657be8;font-weight: 400;">Protein Coding Counts</span> : <%=proteincodingcounts%> / <span style="color:#F9B189;font-weight: 400;">LncRNA Counts</span> : <%=lncrnacounts%></div>
		</div></a>
	</div>
	
	<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Platform :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><%=platform%></div>
		</div>
	</div>

	<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Treatment & Pri/Meta:</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><%=treatment%>  ( <%=primeta%> )</div>
		</div>
	</div>
	

	
		<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Citation(s) :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><%=publication%> ( <%=pubmedid%> )</div>
		</div>
	</div>

    </div>
      
  </div>
</section>

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
</body>
</html>