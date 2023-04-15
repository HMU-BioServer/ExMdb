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

<title>ExMdb : RNA cellular localization</title>



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

div.sp_clear_btn {
    font-size: 5px;
}

div.sp_clear_btn i {
    font-size: 22px;
}

.sp_container {
    width: 180px !important;
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
      <h3>RNA cellular localization</h3>
      <ul>
        <li><a href="index.html">Home</a></li>
        <li><i class='bx bx-chevrons-right'></i></li>
        <li>RNA cellular localization</li>
      </ul>
    </div>
  </div>
  <div class="inner-shape"><img src="assets/images/shape/inner-shape.png" alt="Images"></div>
</div>
<div class="terms-conditions-area" style="padding:40px 0px">
<div class="container wow fadeInUp animated" style="max-width: 1440px;">

		<div class="contact-header">
			<h1>RNA cellular localization</h1>
		</div>		
		<div class="contact-wrapper">
		<div class="row">
		<p class="sub_title_function"><img src="img/iconfont/location.png" style="margin-top: -7px;margin-right: 7px;">ExMdb-celllocation</p>
	    <p class="sub_text_function">Find out the RNA sub-cellular locations for a cell. The identified locations text are shown with white frame.</p>
		<hr>
			<div class="col-lg-2 col-md-2" style="border-right:2px solid #c3c3c3">

				<form id="gene" name="gene" method="post" action="cellmap.jsp" target="exmdb_cell">
				
							<nav class="menu menu--adsila">
								<a class="menu__item" href="#">
									<span class="menu__item-name">Enter key words</span>
								</a>
								<input type="text" name=searchname id="searchname" class="form-control" value="MALAT1" placeholder=" keywords.." style="font-size:20px;">
							</nav>
							<hr>
							<button class="default-btn btn-bg-two border-radius-50" type="submit" style="padding:5px;height:100%;width:90%;font-size:15px" onclick="show_qi()"><span>Submit</span></button>
				</form>
				<nav class="menu menu--adsila">
						<a class="menu__item" href="#">
							<span class="menu__item-name">Workflow :</span>
						</a>
				</nav>
			
				<img src="img/cell/cell_left.png"  style="background: #8eabff40;padding: 10px;border-radius: 10px;">
				
			</div>
			<div class="col-lg-10 col-md-10">
							<iframe frameborder=0 width="100%" height="700px" 
							name="exmdb_cell" id="exmdb_cell"
							src="cellmap.jsp"
							scrolling="no" onload="load()"></iframe>
			</div>
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

<script type="text/javascript" src="select/demo2.js" ></script>
<script type="text/javascript" src="select/jquery.min.js" ></script>
<script type="text/javascript" src="select/selectpage.min.js" ></script>
<script type="text/javascript" src="select/data/cell_out.js" ></script>
<script type="text/javascript">
	$(function(){
		$('#searchname').selectPage({
			showField : 'name',
			keyField : 'id',
			data : tag_data_cell,
			lang: 'en',
			noResultClean : true
		});
		
	
		SyntaxHighlighter.all();
	});
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