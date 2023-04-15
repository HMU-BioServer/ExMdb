<%@ page language="java" import="java.io.*,java.sql.*,java.util.*,java.text.DecimalFormat,wp.base.*" contentType="text/html; charset=UTF-8" pageEncoding="UTF-8"%>
<!DOCTYPE html>
<html>
<head>
<meta charset="ISO-8859-1">
<title>Insert title here</title>

<link rel="stylesheet" href="assets/css/style.css">
<link rel="stylesheet" href="assets/css/bootstrap.min.css">
<link rel="stylesheet" type="text/css" href="tanchu/css/dialog.css" />
<link rel="stylesheet" type="text/css" href="tanchu/css/dialog-donna.css" />
<script src="tanchu/js/modernizr.custom.js"></script>


<style>
body{
font-family:'Montserrat', sans-serif
}
</style>

</head>

<%
String flora="1=1";
String disease="1=1";
String gene="1=1";
String layout_echart="force";

String graph_length="105";
String node_size="20";
String egde_width="2";
String change="0.5";



if(request.getParameter("flora1")!=null&&!request.getParameter("flora1").equals("")){
	flora = request.getParameter("flora1");
}
if(request.getParameter("disease1")!=null&&!request.getParameter("disease1").equals("")){
	disease = request.getParameter("disease1");
}
if(request.getParameter("gene1")!=null&&!request.getParameter("gene1").equals("")){
	gene = request.getParameter("gene1");
}
if(request.getParameter("layout_echart")!=null&&!request.getParameter("layout_echart").equals("")){
	layout_echart = request.getParameter("layout_echart");
}
if(request.getParameter("graph_length")!=null&&!request.getParameter("graph_length").equals("")){
	graph_length = request.getParameter("graph_length");
}
if(request.getParameter("node_size")!=null&&!request.getParameter("node_size").equals("")){
	node_size = request.getParameter("node_size");
}

if(request.getParameter("egde_width")!=null&&!request.getParameter("egde_width").equals("")){
	egde_width = request.getParameter("egde_width");
}

if(request.getParameter("change")!=null&&!request.getParameter("change").equals("")){
	change = request.getParameter("change");
}


String networt_flora=flora;
String network_disease=disease;
String network_gene=gene;

flora = flora.replaceAll(",","','");
disease = disease.replaceAll(",","','");
gene = gene.replaceAll(",","','");

flora = "'"+flora+"'";
disease = "'"+disease+"'";
gene = "'"+gene+"'";

String sql = "select * from exmdb_net where node1 IN ("+flora+","+disease+","+gene+") or node2 IN ("+flora+","+disease+","+gene+")";
String data = dbhello.exmdb_netmaker(sql);

String search_table = dbhello.exmdb_net_table(sql);
dbhello.exmdb_degree_table(sql);
%>



<body>
			<h4 class="Top">Network:</h4>
			<div class="search">
				<p style="padding-bottom: 0px;">
					<font style="font-weight: bold; font-style: italic; color: #070b3b">Your interested Genes</font> Network in EXMDB:
				</p>
				<a data-dialog="somedialog" class="trigger trigger-iframe default-btn btn-bg-two border-radius-50" style="z-index:0;cursor:pointer;padding:7px 20px;text-decoration: none;"> 
					<span>Detailed sub-cellular table of </span> <span style="color:#f5d043;font-weight: 700;text-transform:uppercase">Network</span>
				</a>
				
				<a data-dialog1="somedialog" class="trigger trigger-iframe default-btn btn-bg-two border-radius-50" style="z-index:0;cursor:pointer;padding:7px 20px;text-decoration: none;"> 
					<span>Degree </span> <span style="color:#f5d043;font-weight: 700;text-transform:uppercase">Network</span>
				</a>
				
			</div>
			
				<div id="somedialog" class="dialog" style="z-index: 99;">
					<div class="dialog__overlay"></div>
					<div class="dialog__content" style="height: 700px;border-radius: 10px;">
						<iframe frameborder=0 style="width:100%;height:100%" src="exmdb_local_table.jsp?flora=<%=networt_flora%>&disease=<%=network_disease%>&gene=<%=network_gene%>"></iframe>
						<button class="action default-btn btn-bg-one border-radius-50 ml-20" style="padding:5px 15px;font-family: livvic,sans-serif;" data-dialog-close>Close</button>
					</div>
				</div>
				
				<div id="somedialog_degree" class="dialog" style="z-index: 99;">
					<div class="dialog__overlay"></div>
					<div class="dialog__content" style="height: 700px;border-radius: 10px;">
						<iframe frameborder=0 style="width:100%;height:100%" src="exmdb_degree_table.jsp?flora=<%=networt_flora%>&disease=<%=network_disease%>&gene=<%=network_gene%>"></iframe>
						<button class="action default-btn btn-bg-one border-radius-50 ml-20" style="padding:5px 15px;font-family: livvic,sans-serif;" data-dialog-close>Close</button>
					</div>
				</div>

<script src="tanchu/js/classie.js"></script>
<script src="tanchu/js/dialogFx.js"></script>
<script>
			(function() {

				var dlgtrigger = document.querySelector( '[data-dialog]' ),
					somedialog = document.getElementById( dlgtrigger.getAttribute( 'data-dialog' ) ),
					dlg = new DialogFx( somedialog );

				dlgtrigger.addEventListener( 'click', dlg.toggle.bind(dlg) );

			})();
</script>

<script>
			(function() {

				var dlgtrigger = document.querySelector( '[data-dialog1]' ),
					somedialog = document.getElementById( dlgtrigger.getAttribute( 'data-dialog' ) ),
					dlg = new DialogFx( somedialog_degree );

				dlgtrigger.addEventListener( 'click', dlg.toggle.bind(dlg) );

			})();
</script>

<div id="container_network" style="min-width:400px;height:700px;"></div>





<script type="text/javascript" language="javascript" src="echart/echarts.min.js"></script>

<script>
<%=data%>
</script>

<script type="text/javascript">
        var myChart = echarts.init(document.getElementById('container_network'));
        <%//loading%>
        myChart.showLoading({
            text: 'Loading...',
            textStyle: { fontSize : 30 , color: '#444' },
            effectOption: {backgroundColor: 'rgba(0, 0, 0, 0)'}
        })

        option = {
        	    title: {
        	        text: '',
        	        x: 'center'
        	    },
        	    //鈥斺€  鎮诞妗  鈥斺€ 
        	    tooltip: {
        	    	textStyle:{fontFamily: 'Montserrat',fontSize:15},
        	        formatter: function(x) {
        	        	if(x.color){
        	        		return (x.marker+"Gene Name: "+x.data.name+"<br>"+x.marker+"Gene Type: "+x.data.category); //璁剧疆鎻愮ず妗嗙殑鍐呭鍜屾牸寮  鑺傜偣鍜岃竟閮芥樉绀簄ame灞炴€ 
        	        	}
        	            
        	        },
        	    }, 
        	    legend: [{
        	        orient: 'vertical',
        	        x: 'left',
        	        y: '50px',
        	        itemWidth: 25,
        	        itemHeight: 25,
        	        textStyle:{fontFamily: 'Montserrat',fontSize:15},
        	        data: [ //鑺傜偣鏁版嵁			
        	        	{name:'lincRNA',icon: 'circle'},
        	        	{name:'processed_transcript',icon: 'circle'},
        	        	{name:'antisense',icon: 'circle'},
        	        	{name:'TEC',icon: 'circle'},
        	        	{name:'sense_overlapping',icon: 'circle'},
        	        	{name:'macro_lncRNA',icon: 'circle'},
        	        	{name:'sense_intronic',icon: 'circle'},
        	        	{name:'3prime_overlapping_ncRNA',icon: 'circle'},
        	        	{name:'protein_coding',icon: 'circle'},
        	        	{name:'transcribed_unprocessed_pseudogene',icon: 'circle'},
        	        	{name:'unprocessed_pseudogene',icon: 'circle'},
        	        	{name:'transcribed_processed_pseudogene',icon: 'circle'},
        	        	{name:'transcribed_unitary_pseudogene',icon: 'circle'},
        	        	{name:'processed_pseudogene',icon: 'circle'},
        	        	{name:'pseudogene',icon: 'circle'},
        	        	{name:'unitary_pseudogene',icon: 'circle'},
        	        	{name:'miRNA',icon: 'circle'},
        	        ],
        	        show:true,
        	    }, ],

                toolbox:{
                  show:true,
                  feature:{
                      // 数据视图
                      dataView:{
                          show:true
                      },
                      // 保存图片
                      saveAsImage:{
                          show:true
                      },
                      //动态类型切换
                     // magicType:{
                      //}
                  }
                }, 
        	   
        	    animationDurationUpdate: 30000,
        	    animationEasingUpdate: 'quinticInOut',
        	    
        	    series: [{
        	        type: 'graph',
        	        layout: '<%=layout_echart%>',
        	        circular: {
        	          rotateLabel: true
        	        },
        	        color: ["#FDF5B0","#F2EDB3","#E8E6B7","#DEDFBB","#D4D8BF","#CAD1C3","#C0CAC7","#B6C2CB","#ACBBCF","#A2B4D3","#98ADD7","#8DA6DB","#839FDF","#7998E3","#6F90E7","#6589EB","#5B82EF","#517BF3","#4774F7","#3D6DFB","#3366FF"],
        	        symbolSize: <%=node_size%>,//鍥惧舰鐨勫ぇ灏忥紙绀轰緥涓殑鍦嗙殑澶у皬锛 
        	        roam: true, //榧犳爣缂╂斁鍙婂钩绉 
        	        draggable: true, //鑳藉惁榧犳爣鎷栧姩
        	        focusNodeAdjacency: true, //鏄惁鍦ㄩ紶鏍囩Щ鍒拌妭鐐逛笂鐨勬椂鍊欑獊鍑烘樉绀鸿妭鐐广€佽妭鐐圭殑杈瑰拰閭绘帴鑺傜偣
        	        cursor : 'pointer',
            	    itemStyle: {
          	          normal: {
          	           borderWidth: 1,
          	           borderColor: '#ffffff'
          	           }},
        	        label: {
        	            normal: {
        	                show: true, //鎺у埗闈為珮浜椂鑺傜偣鍚嶇О鏄惁鏄剧ず
        	                position: '',
        	                fontSize: 8,
        	                color: 'black',
        	                fontFamily: 'Montserrat'
        	            },
        	            emphasis: {
        	            	focus: 'adjacency',
        	                lineStyle: {
        	                    width: 10
        	                  },
        	                show: true, //鎺у埗闈為珮浜椂鑺傜偣鍚嶇О鏄惁鏄剧ず
        	                position: 'right',
        	                fontSize: 15,
        	                color: 'black'
        	            },
        	        },
        	        force: {
        	            x: 'center',
        	            y: 'center',
        	            repulsion : 50,
        	            gravity : <%=change%>,
        	            edgeLength : <%=graph_length%>
        	         /* gravity : 0.2,
        	            repulsion: 300,
        	            edgeLength: 100 */
        	            
        	        },
        	        lineStyle:{
        	        	color: 'source',
        	            opacity:0.3,
        	            width:<%=egde_width%>,
        	            type:"solid",
        	            curveness:0.3,
        	            color:"#ccccff",

 					},
        	          emphasis: 
        	          {
        	           lineStyle: {width:10},
        	           itemStyle: {symbolSize: 40}},            
        				data:mydata,
        				links:mylinks,
        	       	  categories: [ 	
        	       		{name:'lincRNA'},
        	       		{name:'processed_transcript'},
        	       		{name:'antisense'},
        	       		{name:'TEC'},
        	       		{name:'sense_overlapping'},
        	       		{name:'macro_lncRNA'},
        	       		{name:'sense_intronic'},
        	       		{name:'3prime_overlapping_ncRNA'},
        	       		{name:'protein_coding'},
        	       		{name:'transcribed_unprocessed_pseudogene'},
        	       		{name:'unprocessed_pseudogene'},
        	       		{name:'transcribed_processed_pseudogene'},
        	       		{name:'transcribed_unitary_pseudogene'},
        	       		{name:'processed_pseudogene'},
        	       		{name:'pseudogene'},
        	       		{name:'unitary_pseudogene'},
        	       		{name:'miRNA'},
        	        ],
        	    }]
        	}; 
        // 使用刚指定的配置项和数据显示图表。
        myChart.setOption(option);
myChart.on('click', function (params) {
	console.log(params)
})
        myChart.hideLoading();
    </script>




</body>
</html>