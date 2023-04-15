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
</style>

<% 




String searchname = "";
if(request.getParameter("searchname")!=null&&!request.getParameter("searchname").equals("")){	
	searchname = request.getParameter("searchname");
}

String sql = "";   
sql = "select * from data_info where i='"+searchname+"'";

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

String down_count = dbhello.exmdb_info_down(searchname);
String up_count = dbhello.exmdb_info_up(searchname);
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
    <div class="section-title">
      <h2 style="margin-left: 0;max-width: 100%;"><%=gseid%>-<span class="in_span"><%=organ%></span> <span style="display:contents;font-size:21px">| Regulatory Genes - Up:[<span class="in_span" style="color:#ff8080"><%=up_count%></span>]/Down:[<span class="in_span" style="color:#009933"><%=down_count%></span>]</span><a style="float:right;font-size: 15px;padding: 10px 15px;" href="exmdb_volcano_out.jsp?GSEID=<%=searchname%>" target="_blank" class="default-btn btn-bg-two border-radius-50">Get Detail <i style="display: contents;font-size: 15px;" class="fa fa-send"></i></a></h2>
      <p class="margin-auto" style="margin-left: 0;max-width: 100%"><span class="sp-color2">Data Processing :</span><%=dataprocessing%></p>
      <p class="margin-auto" style="margin-left: 0;max-width: 100%"><span class="sp-color2">Data Summary :</span><%=summary%></p>
      <p class="margin-auto" style="margin-left: 0;max-width: 100%"><span class="sp-color2">Data Title :</span><%=title%></p>
    <hr>
    </div>  
<div class="contact-wrapper module_box wow fadeInUp animated" style="width:100%">
	<div class="table_contact">
		<div class="row" style="height:100%">
			<div class="col-xl-3 col-lg-3 col-md-3 col-sm-3 left_box">Percentage of samples :</div>
			<div class="col-xl-9 col-lg-9 col-md-9 col-sm-9 right_box"><span style="color: #657be8;font-weight: 400;">Cancer</span> : <%=cancercount%> / <span style="color:#F9B189;font-weight: 400;">Control</span> : <%=normalcount%></div>
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