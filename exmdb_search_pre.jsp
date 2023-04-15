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
<title>Techex - Technology Services HTML Template</title>
</head>

<% 
String searchname = "HOTAIR";  



if(request.getParameter("searchname")!=null&&!request.getParameter("searchname").equals("")){
	searchname = request.getParameter("searchname");
}



%>


<body>
<div class="preloader">
  <div class="d-table">
    <div class="d-table-cell">
      <div class="spinner"></div>
    </div>
  </div>
</div>
<header class="top-header top-header-bg">
  <div class="container">
    <div class="row align-items-center">
      <div class="col-lg-7 col-md-6">
        <div class="top-head-left">
          <div class="top-contact">
            <h3>Support By:<a href="">+1(212) 255-5511</a></h3>
          </div>
        </div>
      </div>
      <div class="col-lg-5 col-md-6">
        <div class="top-header-right">
          <div class="top-header-social">
            <ul>
              <li><a href="https://www.facebook.com/" target="_blank"><i class='bx bxl-facebook'></i></a></li>
              <li><a href="https://twitter.com/?lang=en" target="_blank"><i class='bx bxl-twitter'></i></a></li>
              <li><a href="https://www.linkedin.com/" target="_blank"><i class='bx bxl-linkedin-square'></i></a></li>
              <li><a href="https://www.instagram.com/" target="_blank"><i class='bx bxl-instagram'></i></a></li>
            </ul>
          </div>
          <div class="language-list">
            <select class="language-list-item">
              <option>English</option>
              <option>Ø§ÙØ¹Ø±Ø¨ÙÙØ©</option>
              <option>Deutsch</option>
              <option>PortuguÃªs</option>
              <option>ç®ä½ä¸­æ</option>
            </select>
          </div>
        </div>
      </div>
    </div>
  </div>
</header>
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
      <h3>Terms & Conditions</h3>
      <ul>
        <li><a href="index.html">Home</a></li>
        <li><i class='bx bx-chevrons-right'></i></li>
        <li>Terms & Conditions</li>
      </ul>
    </div>
  </div>
  <div class="inner-shape"><img src="assets/images/shape/inner-shape.png" alt="Images"></div>
</div>
<div class="terms-conditions-area pt-100 pb-70">
<div class="row">

</div>
<div class="" style="width:60%;margin:auto" >


<h1>Body Map</h1>
						<iframe frameborder=0 width="100%" height="850px"
							name="exmdb_bodymap" id="exmdb_bodymap"
							src="bodymap.jsp?searchname=<%=searchname%>"
							scrolling="no" ></iframe>

<h1>Cell Map</h1>
						<iframe frameborder=0 width="100%" height="700px"
							name="exmdb_cellmap" id="exmdb_cellmap"
							src="cellmap.jsp?searchname=<%=searchname%>"
							scrolling="no" ></iframe>
							
<h1>Exp data</h1>
						<iframe frameborder=0 width="100%" height="700px"
							name="exmdb_exp" id="exmdb_exp"
							src="exmdb_exp_data_table.jsp?searchname=<%=searchname%>"
							scrolling="no" ></iframe>
<h1>Biomarker</h1>
						<iframe frameborder=0 width="100%" height="700px"
							name="exmdb_bio" id="exmdb_bio"
							src="exmdb_biomarker_data_table.jsp?searchname=<%=searchname%>"
							scrolling="no" ></iframe>

</div>

</div>
<footer class="footer-area footer-bg">
  <div class="container">
    <div class="footer-top pt-100 pb-70">
      <div class="row">
        <div class="col-lg-4 col-sm-6">
          <div class="footer-widget">
            <div class="footer-logo"><a href="index.html"><img src="assets/images/logos/footer-logo.png" alt="Images"></a></div>
            <p>Proin gravida nibh vel velit auctor aliquet. Aenean sollicitudin,lorem quis bibendum auct.Aenean,lorem quis bibendum auct. Aenean sollicitudin lorem. </p>
            <div class="footer-call-content">
              <h3>Talk to Our Support</h3>
              <span><a href="">+1 002-123-4567</a></span><i class='bx bx-headphone'></i></div>
          </div>
        </div>
        <div class="col-lg-2 col-sm-6">
          <div class="footer-widget pl-2">
            <h3>Services</h3>
            <ul class="footer-list">
              <li><a href="service-details.html" target="_blank"><i class='bx bx-chevron-right'></i>IT Consultancy </a></li>
              <li><a href="service-details.html" target="_blank"><i class='bx bx-chevron-right'></i>Business Solution </a></li>
              <li><a href="service-details.html" target="_blank"><i class='bx bx-chevron-right'></i>Digital Services </a></li>
              <li><a href="compare.html" target="_blank"><i class='bx bx-chevron-right'></i>Business Reform </a></li>
              <li><a href="service-details.html" target="_blank"><i class='bx bx-chevron-right'></i>Web Development </a></li>
              <li><a href="service-details.html" target="_blank"><i class='bx bx-chevron-right'></i>Cloud Computing </a></li>
              <li><a href="service-details.html" target="_blank"><i class='bx bx-chevron-right'></i>Data Analysis </a></li>
            </ul>
          </div>
        </div>
        <div class="col-lg-3 col-sm-6">
          <div class="footer-widget pl-5">
            <h3>Our Blog</h3>
            <ul class="footer-blog">
              <li><a href="blog-details.html"><img src="assets/images/blog/blog-img-footer.jpg" alt="Images"></a>
                <div class="content">
                  <h3><a href="blog-details.html">Product Idea Solution For New Generation</a></h3>
                  <span>04 Dec 2020</span></div>
              </li>
              <li><a href="blog-details.html"><img src="assets/images/blog/blog-img-footer2.jpg" alt="Images"></a>
                <div class="content">
                  <h3><a href="blog-details.html">New Device Invention for Digital Platform</a></h3>
                  <span>07 Dec 2020</span></div>
              </li>
              <li><a href="blog-details.html"><img src="assets/images/blog/blog-img-footer3.jpg" alt="Images"></a>
                <div class="content">
                  <h3><a href="blog-details.html">Business Strategy Make His Goal Acheive</a></h3>
                  <span>10 Dec 2020</span></div>
              </li>
            </ul>
          </div>
        </div>
        <div class="col-lg-3 col-sm-6">
          <div class="footer-widget">
            <h3>Newsletter</h3>
            <p>Lorem ipsum dolor sit amet,consectetur adipiscing elit. Vestibulum finibus molestie molestie. Phasellus ac rutrum massa,et volutpat nisl. Fusce ultrices suscipit nisl.</p>
            <div class="newsletter-area">
              <form class="newsletter-form" data-toggle="validator" method="POST">
                <input type="email" class="form-control" placeholder="Enter Your Email" name="EMAIL" required autocomplete="off">
                <button class="subscribe-btn" type="submit"><i class='bx bx-paper-plane'></i></button>
                <div id="validator-newsletter" class="form-result"></div>
              </form>
            </div>
          </div>
        </div>
      </div>
    </div>
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
</body>
</html>