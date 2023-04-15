<%@ page language="java" import="java.sql.*,java.util.*,wp.*" contentType="text/html; charset=UTF-8"
    pageEncoding="UTF-8"%>
<!DOCTYPE html>
<html lang="zxx">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta http-equiv="X-UA-Compatible" content="IE=edge" />
    <meta http-equiv="X-UA-Compatible" content="ie=edge">
    <title>CellTracer: Help</title>
	
	<!-- small icon -->
	<link rel="shortcut icon" href="my_images/icon_logo.png">

      <!-- Stylesheet -->
    <link rel="stylesheet" href="assets/css/vendor.css">
    <link rel="stylesheet" href="assets/css/nice-select.css">
    <link rel="stylesheet" href="assets/css/magnific-popup.css">
    <link rel="stylesheet" href="assets/css/fancybox.min.css">
    <link rel="stylesheet" href="assets/css/style.css">
    <link rel="stylesheet" href="assets/css/responsive.css">

    
     <!-- animate -->
    <link rel="stylesheet" href="Animate/animate.min.css">
    <script src="Animate/jquery-1.8.3.min.js"></script>
	<script src="http://www.jq22.com/jquery/jquery-migrate-1.2.1.min.js"></script>

<style>

.card-body{
text-align:justify;
}
</style>
</head>
<body>

    <!-- preloader area start -->
    <div class="preloader" id="preloader">
        <div class="preloader-inner">
            <div class="spinner">
                <div class="dot1"></div>
                <div class="dot2"></div>
            </div>
        </div>
    </div>
    <!-- preloader area end -->
   

   
        <!-- navbar start -->
    

      <div class="navbar-area navbar-area-2">
        <nav class="navbar navbar-expand-lg">
            <div class="container nav-container">
                
                <div class="logo">
                    <a class="main-logo" href="CellTracer_index.jsp"><img src="my_images/logo.png" alt="img" style="width:250px" ></a>
                </div>
                <div class="nav-right-part nav-right-part-mobile">
                    <a class="btn btn-base" href="#">Login</a>
                </div>
                <div class="collapse navbar-collapse" id="dkt_main_menu">
                    <ul class="navbar-nav menu-open">
                        <li >
                            <a href="CellTracer_index.jsp">Home</a>                            
                        </li>
                      <li class="menu-item-has-children current-menu-item">
                            <a href="#">Search</a>    
                             <ul class="sub-menu" >
                              <li ><a href="CellTracer_Search_gene.jsp">Search by gene</a></li>
                                 <li ><a href="CellTracer_Search_disease.jsp">Search by disease</a></li>
                                <li ><a href="CellTracer_Search_organization.jsp">Search by organ</a></li>
                                <li ><a href="CellTracer_Search_species.jsp">Search by species</a></li>  
                                <li ><a href="CellTracer_Search_blast.jsp">Search by sequence</a></li>                               
                                 <li ><a href="CellTracer_Search_customized.jsp">Search by  customize settings</a></li>                                     
                            </ul>                                                                             
                        </li>
                        <li class="menu-item-has-children current-menu-item">
                            <a href="#">Tools</a>
                            <ul class="sub-menu" style="left: -270%;">
                              <li><a style="background: #ffb980;color: white;text-align: center;font-size:18px">Comprehensive tools</a></li>
                                 <li class="colorControl"><a href="CellTracer_ComprehensiveTool-1.jsp">CellTracer- geneCellCluster</a></li>
                                <li class="colorControl"><a href="CellTracer_ComprehensiveTool-2.jsp">CellTracer- geneCellTraject</a></li>
                                 <li class="colorControl"><a href="CellTracer_ComprehensiveTool-3.jsp">CellTracer- geneStateTraject</a></li>  
                                 <li class="colorControl"><a href="CellTracer_ComprehensiveTool-4.jsp">CellTracer- geneFuncTraject</a></li> 
                                    <li class="colorControl"><a href="CellTracer_3D.jsp">CellTracer- Multi-Omics-3D</a></li>                                                                          
                            </ul>
                            
                              <ul class="sub-menu"  style="left: 130%;">
                                <li class="colorControl1"><a style="    background: #5ab1ef;color: white;text-align: center;font-size:18px">Mini analysis tools</a></li>
                                 <li class="colorControl1"><a href="CellTracer_Cluster.jsp">CellTracer- CellCluster</a></li>
                                <li class="colorControl1"><a href="CellTracer_Trajectory.jsp">CellTracer- CellTraject</a></li>
                                 <li class="colorControl1"> <a href="CellTracer_Timeline.jsp">CellTracer- CellTimeLine</a></li>                                
                                <li class="colorControl1"><a href="CellTracer_GeneExp.jsp">CellTracer- GeneExpression</a></li>                         
                                 <li class="colorControl1"><a href=CellTracer_Survival.jsp>CellTracer- GeneSurvival</a></li>                              
                                <li class="colorControl"><a href="CellTracer_CellCor.jsp">CellTracer- GeneStateInterplay</a></li>
                                <li class="colorControl1"><a href="CellTracer_PathwayCor.jsp">CellTracer- GeneFuncInterplay</a></li>                                                
                            </ul>
                            
                        </li>
                         <li >
                            <a href="CellTracer_Download.jsp">Download</a>                            
                        </li>
                        <li >
                            <a href="CellTracer_Statistic.jsp">Statistic</a>                            
                        </li>
                       
                         <li><a href="CellTracer_Help.jsp">Help</a>
                        </li>
                    </ul>
                </div>
               
            </div>
        </nav>
    </div>
    <!-- navbar end -->
    
    <!-- breadcrumb start -->
    <div class="breadcrumb-area" style="background: url(my_images/banner2.jpg) no-repeat center -100px;background-size: cover;">
        <div class="container">
            <div class="row">
                <div class="col-lg-12">
                    <div class="breadcrumb-inner">
                        <div class="section-title text-center">
                            <h2 class="page-title">FAQ</h2>
                            <ul class="page-list">
                                <li><a href="CellTracer_index.jsp">Home</a></li>
                                <li>FAQ</li>
                            </ul>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>
    <!-- breadcrumb end -->

        <!-- faq-area start -->
    <section class="faq-page-area pd-top-100 pd-bottom-100">
        <div class="container">
            <div class="row justify-content-center">
                <div class="col-lg-12">
                    <div class="faq-accordion accordion" id="accordionExample">
                    
                     <div class="card">
                            <div class="card-header">
                                <h2>
                                    <button class="btn collapsed" type="button" data-toggle="collapse" data-target="#collapseZero">
                                  <span style="color: #ffb980;font-weight: bold;">Cell</span><span style="color: #5ab1ef;font-weight: bold;">Tracer</span> Introduction
                                    <span class="collapse-icon"></span>
                                    </button>
                                </h2>
                            </div>
                            <div id="collapseZero" class="collapse show" data-parent="#accordionExample">
                                <div class="card-body">
								<p style="text-align: center;font-weight: bold;display: block;">CellTracer: A comprehensive database to dissect causative multi-omics interplay contributing to cellular development trajectory</p>		
<br>

<p>The current release of CellTracer database contains <b>194,1552</b> cells from <b>222</b> single-cell datasets of <b>42</b> diseases and <b>80</b> organs. A number of <b>117</b> cell types were identified from different cellular microenvironment such as clinical treatments of Chemotherapy, Immunotherapy or other Targeted therapy methods.</p>
<br>
<p>Based on CellTracer, users can explore the causative multi-omics interplay between more than <b>50,000</b> genes (including coding-genes, lncRNAs, pseodugenes and etc.) and more than <b>10,000</b> biological contexts (including Gene Ontologies, Pathways, Cancer Hallmarks, Cellular States and etc. ) contributing to cellular development trajectory.</p>

                                </div>
                            </div>
                        </div>  
                    
                    
                        <div class="card">
                            <div class="card-header">
                                <h2>
                                    <button class="btn collapsed" type="button" data-toggle="collapse" data-target="#collapseZeroo">
                                     <span style="color: #ffb980;font-weight: bold;">Cell</span><span style="color: #5ab1ef;font-weight: bold;">Tracer</span>'s Home
                                    <span class="collapse-icon"></span>
                                    </button>
                                </h2>
                            </div>
                            <div id="collapseZeroo" class="collapse" data-parent="#accordionExample">
                                <div class="card-body">
								<p style="text-align: center;font-weight: bold;display: block;">CellTracer provides quick search and powerful analysis tools and high-throughput single-cell data datasets</p><br>

<p><b>1.</b> Click this search button to start a quick search for gene/organ/disease etc.</p>

<p><b>2.</b> Click this button to get an introduction to CellTracer.</p>

<p><b>3.</b> CellTracer's navigation bar.</p>

<p><b>4.</b> Click this button to get our contact information and our recent works.</p>

<p><b>5.</b> All of CellTracer's powerful analysis tools.</p>

<p><b>6.</b> CellTracer's statistical data results.</p>

	<img style="  margin-top: 50px;  height:480px" src="help_img/h2-1.png">
	
	
	<img style="width:1040px" src="help_img/h2-2.png">

<p style="text-align:center">
			
			<span style="display: inline-table;width: 50%;">Figure 2-1</span>
			</p>		
                                </div>
                            </div>
                        </div>  
                    
                    
                        <div class="card">
                            <div class="card-header">
                                <h2>
                                    <button class="btn collapsed" type="button" data-toggle="collapse" data-target="#collapseOne">
                                    Quick Search in  <span style="color: #ffb980;font-weight: bold;">Cell</span><span style="color: #5ab1ef;font-weight: bold;">Tracer</span>
                                    <span class="collapse-icon"></span>
                                    </button>
                                </h2>
                            </div>
                            <div id="collapseOne" class="collapse" data-parent="#accordionExample">
                                <div class="card-body">
								<p style="text-align: center;font-weight: bold;display: block;"><span style="color:#F28705">CellTracer</span> provides a user-friendly search and result interface.</p><br>

<p><b>1.</b> You can search for multiple data types in CellTracer, and CellTracer will return the results of the relevant dataset in the background based on your input.</p>

<p><b>2.</b> Enter what you want to find in the search box.</p>
<img style="  margin-top: 50px;  height:250px" src="help_img/h3-1.png">
<p style="text-align:center">
			
			<span style="display: inline-table;width: 50%;">Figure 3-1</span>
			</p>	
                                </div>
                            </div>
                        </div>                      
                        <div class="card">
                            <div class="card-header">
                                <h2>
                                    <button class="btn collapsed" type="button" data-toggle="collapse" data-target="#collapseTwo">
                                    
                                    Find your interested datasets in  <span style="color: #ffb980;font-weight: bold;">Cell</span><span style="color: #5ab1ef;font-weight: bold;">Tracer</span>
                                    <span class="collapse-icon"></span>
                                    </button>
                                </h2>
                            </div>
                            <div id="collapseTwo" class="collapse" data-parent="#accordionExample">
								<div class="card-body">
								
<p style="text-align: center;font-weight: bold;display: block;">CellTracer provides a variety of search methods including gene names, disease names, tissues, species, sequences and customized searches (For example: customized search).</p><br>
<p><b>1.</b> Select the data type you are interested in to locate the data set in CellTracer.</p>
<p><b>2.</b> Result display.</p>
								<img style="  margin-top: 50px;  height:250px" src="help_img/h4-1.png">
									<img style=" margin-top: 40px;   height:250px;width: 680px;" src="help_img/h4-2.png">
			
			<p style="display: flex;">
			<span style="display: inline-table;width: 50%;padding-left: 13%;">Figure 4-1</span>
			<span style="display: inline-table;width: 50%;padding-left: 15%;">Figure 4-2</span>
			</p>						
									
								</div>
							</div>
                        </div>
                       
                       
                        <div class="card">
                            <div class="card-header">
                                <h2>
                                    <button class="btn collapsed" type="button" data-toggle="collapse" data-target="#collapseThree">
                                   Basic information for Dataset / Cell in <span style="color: #ffb980;font-weight: bold;">Cell</span><span style="color: #5ab1ef;font-weight: bold;">Tracer</span>
                                    <span class="collapse-icon"></span>
                                    </button>
                                </h2>
                            </div>
                            <div id="collapseThree" class="collapse" data-parent="#accordionExample">
                                <div class="card-body">
<p style="text-align: center;font-weight: bold;display: block;">For each dataset or cell unit, <span style="color:#F28705">CellTracer</span> provides basic information, which is displayed in Figure 5-1 and Figure 5-2.</p><br>
<p><b>1.</b> The basic information of the dataset event you searched.</p>
<p><b>2.</b> You can click on the analysis results of different presentations on the right to further explore the data set that interests you.</p>

<p><b>3.</b> The basic information of the cell event you searched.</p>

<img style="  margin-top: 50px;  " src="help_img/h5-1.png">

<p style="text-align:center">
			
			<span style="display: inline-table;width: 50%;">Figure 5-1</span>
			</p>	
									<img style=" margin-top: 40px; " src="help_img/h5-2.png">
			
					<p style="text-align:center">
			
			<span style="display: inline-table;width: 50%;">Figure 5-2</span>
			</p>	


</div>
                            </div>
                        </div>
                        <div class="card">
                            <div class="card-header">
                                <h2>
                                    <button class="btn collapsed" type="button" data-toggle="collapse" data-target="#collapseFour">
                                  5 Comprehensive analysis tools in  <span style="color: #ffb980;font-weight: bold;">Cell</span><span style="color: #5ab1ef;font-weight: bold;">Tracer</span>
                                    <span class="collapse-icon"></span>
                                    </button>
                                </h2>
                            </div>
                            <div id="collapseFour" class="collapse" data-parent="#accordionExample">
                                <div class="card-body">
<p style="text-align: center;font-weight: bold;display: block;">CellTracer provides 5 comprehensive analysis tools: geneCellCluster, geneCellTraject, geneStateTraject, geneFuncTraject, Multi-Omics-3D.You can use these tools to explore single-cell datasets in depth (For example: geneCellCluster).</p><br>
<p><b>1.</b> CellTracer also allows users to customize chart parameters to get the best results after analyzing single-cell datasets. You can set point size, search for different gene names, cluster methods, different resolutions, and point coordinates (including tSNE and UMAP).</p>
<p><b>2.</b> Visualization of CellTracer analysis results (including presenting the results of different classifications of this dataset,  showing the gene expression values in each cell of the dataset and displaying the average expression value of genes in each cluster etc.).</p>

					<img style=" margin-top: 30px; " src="help_img/h6-1.png">
			
			<p style="text-align:center">
			
			<span style="display: inline-table;width: 50%;">Figure 6-1</span>
			</p>					
</div>
                            </div>
                        </div>
                        <div class="card">
                            <div class="card-header">
                                <h2>
                                    <button class="btn collapsed" type="button" data-toggle="collapse" data-target="#collapseFive">
                                     7 Mini analysis tools in  <span style="color: #ffb980;font-weight: bold;">Cell</span><span style="color: #5ab1ef;font-weight: bold;">Tracer</span>
                                    <span class="collapse-icon"></span>
                                    </button>
                                </h2>
                            </div>
                            <div id="collapseFive" class="collapse" data-parent="#accordionExample">
                                <div class="card-body">
<p style="text-align: center;font-weight: bold;display: block;">CellTracer provides 7 mini analysis tools: CellCluster, CellTraject, CellTimeLine, GeneExpression, GeneSurvival, GeneStateInterplay, GeneFuncInterplay. If you don't want to do in-depth analysis of the dataset, so you can use a fast mini analysis of the single-cell data set to get the desired results (For example: GeneSurvival).</p><br>
<p><b>1.</b> You can select the data and genes you are interested in, and you can set the color and thickness of the lines and how to differentiate between high and low expression groups to get the desired results.</p>
<p><b>2.</b> GeneSurviva's analysis result.</p>


					<img style=" margin-top: 30px;   height:450px;width:50%" src="help_img/h7-1.png">
				<img style=" margin-top: 30px;   height:450px;width:49%" src="help_img/h7-2.png">
			<p style="display: flex;">
			<span style="display: inline-table;padding-left: 22%;">Figure 7-1</span>
			<span style="display: inline-table;padding-left: 45%;">Figure 7-2</span>
			</p>						
</div>
                            </div>
                        </div>
                        <div class="card">
                            <div class="card-header">
                                <h2>
                                    <button class="btn collapsed" type="button" data-toggle="collapse" data-target="#collapseSix">
                                    Download in  <span style="color: #ffb980;font-weight: bold;">Cell</span><span style="color: #5ab1ef;font-weight: bold;">Tracer</span>
                                    <span class="collapse-icon"></span>
                                    </button>
                                </h2>
                            </div>
                            <div id="collapseSix" class="collapse" data-parent="#accordionExample">
                                <div class="card-body">
<p style="text-align: center;font-weight: bold;display: block;">How to download dataset.</p><br>
<p><b>1.</b> We provide R language data analysis results for 222 single cell datasets. You can click Datasize to download the results for secondary utilization of the datasets. We hope CellTracer can fill the gaps in the single cell database.</p>


					<img style=" margin-top: 30px;   height:450px;width:100%" src="help_img/h8-1.png">
			
			<p style="display: flex;">
			<span style="display: inline-table;padding-left: 45%;">Figure 8-1</span>
			</p>						
</div>
                            </div>
                        </div>
                       
                 
                    </div>
                </div>  
            </div>
        </div>          
    </section>
    <!-- faq-area start -->
    <!-- footer area start -->
   <footer class="footer-area ">       
        <div class="container">
            <div class="copyright-area">
                <div class="row" style="    text-align: center;">
                    <div class="col-lg-12 align-self-center">
                        <p>© College of Bioinformatics Science and Technology, Harbin Medical University, Harbin, China.
   <svg style="display: inline;" t="1624811063253" class="icon" viewBox="0 0 1024 1024" version="1.1" xmlns="http://www.w3.org/2000/svg" p-id="18328" width="30" height="30"><path d="M512 512m-512 0a512 512 0 1 0 1024 0 512 512 0 1 0-1024 0Z" fill="#E94B35" p-id="18329"></path><path d="M266.666667 373.333333l-94.037334 65.450667 33.173334-109.653333-91.306667-69.248 114.56-2.346667L266.666667 149.333333l37.610666 108.224 114.56 2.346667-91.306666 69.205333 33.173333 109.653334zM442.154667 186.709333l-29.525334 7.808 18.304-24.448-16.554666-25.664 28.906666 9.856 19.306667-23.68-0.426667 30.549334 28.48 11.029333-29.205333 9.024-1.685333 30.506667zM541.866667 268.693333l-13.013334 27.648-6.933333-29.738666-30.314667-3.84 26.133334-15.786667L512 216.96l23.125333 19.968 26.752-14.72-11.84 28.16 22.272 20.906667zM530.368 381.546667l-27.733333 12.821333 13.802666-27.264-20.778666-22.4 30.186666 4.693333 14.890667-26.666666 4.864 30.144 29.994667 5.930666-27.2 13.930667 3.626666 30.336zM442.154667 464.042667l-29.525334 7.808 18.304-24.448-16.554666-25.664 28.906666 9.856 19.306667-23.68-0.426667 30.549333 28.48 11.029333-29.205333 9.024-1.685333 30.506667z" fill="#F2C500" p-id="18330"></path></svg>
                        
                        
                        </p>
                    </div>
                   
                </div>
            </div>                
        </div>
    </footer>     
    <!-- footer area end -->

    <!-- back to top area start -->
    <div class="back-to-top">
        <span class="back-top"><i class="fa fa-angle-up"></i></span>
    </div>
    <!-- back to top area end -->


    <!-- all plugins here -->
    <script src="assets/js/vendor.js"></script>
    <!-- main js  -->
    <script src="assets/js/main.js"></script>
</body>
</html>