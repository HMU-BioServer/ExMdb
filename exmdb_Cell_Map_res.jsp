<%@ page language="java" import="java.io.*,java.sql.*,java.util.*,java.text.DecimalFormat,wp.base.*" contentType="text/html; charset=UTF-8"
    pageEncoding="UTF-8"%>
<!DOCTYPE html>
<html>
<head>
<meta charset="ISO-8859-1">



<title>CeRNA-Cell-Map</title>
<!-- Main style sheet -->

<link rel="stylesheet" href="assets/css/bootstrap.min.css">

<link rel="stylesheet" href="select/selectpage_cell.css" type="text/css">
<!-- echarts -->
<script type="text/javascript" src="echart/echarts.min.js"></script>
  

</head>
<%


String location="umap";
if(request.getParameter("location")!=null&&!request.getParameter("location").equals("")){	
	location = request.getParameter("location");
}

// System.out.println(location);

String resolution="5";
if(request.getParameter("resolution")!=null&&!request.getParameter("resolution").equals("")){	
	resolution = request.getParameter("resolution");
}

String dot_2_type = "categories:lable";
String color_bar = "'rgba(11,39,81,0.8)','rgb(47, 47, 183)','rgba(97,156,255,0.8)', 'rgba(0,191,196,0.5)', 'rgba(57,182,0,0.7)','rgba(255,255,72,0.9)','rgb(255, 102, 0)','rgba(255,98,188,0.5)','rgba(255,98,188,0.7)','rgba(204,0,0,0.6)','rgba(204,0,0,0.9)','rgb(255, 0, 102)','rgb(255, 0, 0)'";

String gseid="BLCA_GSE145281_aPDL1_count";
if(request.getParameter("gseid")!=null&&!request.getParameter("gseid").equals("")){	
	gseid = request.getParameter("gseid");
}
// String data = dbhello.getCelltSNE("select * from tsne_cell");

String character="sample";

String pseudotime_max_min="_";
String pseudotime_max="0";
String pseudotime_min="0";

if(request.getParameter("character")!=null&&!request.getParameter("character").equals("")){	
	character = request.getParameter("character");
	if (character.equals("pseudotime")){
		pseudotime_max_min = dbhello.exmdb_get_pseudotime_maxmin(gseid);
		String Linshi [] = pseudotime_max_min.split("_");
		pseudotime_max = Linshi[0];
		pseudotime_min = Linshi[1];
		dot_2_type = "";
		color_bar = "'rgb(255, 224, 204)','rgb(230, 90, 0)'";
		
	}
}

// System.out.println("max: "+pseudotime_max);
// System.out.println("min: "+pseudotime_min);

String searchname="MALAT1";
if(request.getParameter("searchname")!=null&&!request.getParameter("searchname").equals("")){	
	searchname = request.getParameter("searchname");
}


String path = request.getRealPath("/").replace("\\","/");   //get root real path 
String immpath  = path+"sc_exp/"+gseid+".txt";
BufferedReader br = new fcon(immpath).getBr();
String str = br.readLine(); // igron 1 line


int cellcount = dbhello.exmdb_get_cell_count(gseid);
double QHT_array[] = new double[cellcount];

// System.out.println("cellcount: " + cellcount);


while(str!=null && str.length()>0) {
	if(str.startsWith(searchname)) {
		String geneidsArray [] = str.split(" ");
		String gene_index [] = geneidsArray[1].split("#");
// 		System.out.println("for_out" + " length = " + gene_index.length);
			for(int i=0;i<gene_index.length;i++){
				String index_exp[] =  gene_index[i].split("_");
// 				System.out.println(Integer.parseInt(index_exp[0]));
// 				System.out.println(gene_index[i]);
				QHT_array[(Integer.parseInt(index_exp[0])-1)] = Double.parseDouble(index_exp[1]);
			}
		break;
	}
	str = br.readLine();
}



String data = dbhello.getCelltSNE_new("select * from tsne_cell_all where dataid = '"+gseid+"'",location,resolution,character,QHT_array);
// String data = dbhello.getCelltSNE("select * from tsne_cell_all");

String location_1 = "";
String label = "";
String max="10";
String min="1";
String strb = dbString.getSelectItem1("select concat(max("+"resolution"+resolution+"_class"+"),',',min("+"resolution"+resolution+"_class"+")) from tsne_cell_all where dataid='"+gseid+"'");


max = strb.split(",")[0];
min = strb.split(",")[1];


// System.out.println("location: "+location+"\n"+"resolution: "+resolution+"\n"+"gseid: "+gseid+"\n"+"character :"+character);

%> 
<body>


<div class="row">
<div class="col-xl-4 col-lg-4 col-md-4 col-sm-4" id="container_left" style="min-width:400px;height:400px;"></div>
<div class="col-xl-4 col-lg-4 col-md-4 col-sm-4" id="container_right" style="min-width:400px;height:400px;"></div>
<div class="col-xl-4 col-lg-4 col-md-4 col-sm-4" id="container_right_right" style="min-width:400px;height:400px;"></div>
</div>


<!-- left --> 

 <script type="text/javascript">
 <%=data%>
 </script>
     
 <script type="text/javascript">
        // 基于准备好的dom，初始化echarts实例
        var myChart = echarts.init(document.getElementById('container_left'));
        <%//loading%>
        myChart.showLoading({
            text: 'Loading...',
            textStyle: { fontSize : 30 , color: '#444' },
            effectOption: {backgroundColor: 'rgba(0, 0, 0, 0)'}
        })
        option = {
                title: {
                    text: 'Cell clusters',
                    textStyle:{fontFamily: 'Montserrat'},
                },
                
                 // 工具箱
              toolbox:{
                show:true,
                feature:{
                    // 数据视图
                    dataView:{
                        show:true
                    },
                    // 还原
                    restore:{
                        show:true
                    },
                    // 区域缩放
                    dataZoom:{
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
              
              
                visualMap: {
                    dimension:2,
                    min: <%=min%>,           //----------------------------
                    max: <%=max%>,         //----------------------------
                    orient: 'vertical',
                    right: 0,
                    textStyle:{fontFamily: 'Montserrat'},
                    top: 'center',
                    text: ['Clusters', ''],
                    calculable: true,
                    inRange: {
                        color: [
                        'rgba(11,39,81,0.8)',
                        'rgba(97,156,255,0.8)', 
                        'rgba(0,191,196,0.5)', 
                        'rgba(57,182,0,0.7)',
                        'rgba(255,255,72,0.9)',
                        'rgba(255,98,188,0.5)',
                        'rgba(255,98,188,0.7)',
                        'rgba(204,0,0,0.6)',
                        'rgba(204,0,0,0.9)',
                        ]
                    },
                },
                
                tooltip: {
                	textStyle:{fontFamily: 'Montserrat'},
                    trigger: 'item',
                    axisPointer: {
                        type: 'cross'
                    },
                    
                 //   formatter:'{c}'
                    
                    formatter: function (params) {
                        //params.value data= [[1,1,'a',5],[1,1,'a',5]]
                        //return 'Cell:'+params.value[2]+'<br> Cluster: '+params.value[3]
                        //data=[{name:'aaa',value:[1,1,3]},{name:'aaa',value:[2,2,2]},],
                        //console.log(params.data.name)
                        //console.log(params.data.value[0])//value:[1,1,3]中第一个数
                        return params.marker+'Cell:'+params.name+'<br>'+params.marker+'Cluster: '+params.value[2]
                        
                    }
                },
                grid: {
                	top:'13%',
                    left: '5%',
                    right: '18%',
                    bottom: '5%',
                    containLabel: true
                  },
                
                xAxis: [{
                	textStyle:{fontFamily: 'Montserrat'},
                name:"<%=location%>1",
                nameLocation:"middle",//middle start end
                nameGap:20,
                type: "value",
               // data: xData,当type为value时,data不起作用，此时使用的时series中的第一维数据
            }],
            
            
               yAxis: [{
            	   textStyle:{fontFamily: 'Montserrat'},
                type: "value",
                name:"<%=location%>2",
                nameLocation:"middle",//middle start end
                nameGap:20,
            }],
            
            
                series: [{
                //    name: 'price-area',
                    type: 'scatter',
                    symbolSize: 2,
                     itemStyle: {
                         normal: {
                             borderWidth: 0.5,
                             borderColor: '#fff',
                         }
                     },
                    data: location_1
                }]
            };

        // 使用刚指定的配置项和数据显示图表。
        myChart.setOption(option);

        <%//loading%>
        myChart.hideLoading();
    </script>
    

<script type="text/javascript">
        // 基于准备好的dom，初始化echarts实例
        var myChart = echarts.init(document.getElementById('container_right'));
        <%//loading%>
        myChart.showLoading({
            text: 'Loading...',
            textStyle: { fontSize : 30 , color: '#444' },
            effectOption: {backgroundColor: 'rgba(0, 0, 0, 0)'}
        })
        option = {
                title: {
                    text: 'Cell type',
                    textStyle:{fontFamily: 'Montserrat'},
                },
                
                 // 工具箱
              toolbox:{
                show:true,
                feature:{
                    // 数据视图
                    dataView:{
                        show:true
                    },
                    // 还原
                    restore:{
                        show:true
                    },
                    // 区域缩放
                    dataZoom:{
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
              
              
                visualMap: {
                    min: <%=pseudotime_min%>,           //----------------------------
                    max: <%=pseudotime_max%>,         //----------------------------
                    orient: 'vertical',
                    right: 0,
                    textStyle:{fontFamily: 'Montserrat'},
                    top: 'center',
                    calculable: true,
                    inRange: {
                        color: [
                        <%=color_bar%>
                        ]
                    },
                    <%=dot_2_type%>
                },
                
                tooltip: {
                	textStyle:{fontFamily: 'Montserrat'},
                    trigger: 'item',
                    axisPointer: {
                        type: 'cross'
                    },
                    
                 //   formatter:'{c}'
                    
                    formatter: function (params) {
                        //params.value data= [[1,1,'a',5],[1,1,'a',5]]
                        //return 'Cell:'+params.value[2]+'<br> Cluster: '+params.value[3]
                        //data=[{name:'aaa',value:[1,1,3]},{name:'aaa',value:[2,2,2]},],
                        //console.log(params.data.name)
                        //console.log(params.data.value[0])//value:[1,1,3]中第一个数
                        return params.marker+'Cell:'+params.name+'<br>'+params.marker+"<%=character%>"+": "+params.value[4]
                        
                    }
                },
                legend: {

              	  type: 'scroll'
                },
                grid: {
                	top:'13%',
                    left: '5%',
                    right: '15%',
                    bottom: '5%',
                    containLabel: true
                  },
                
                xAxis: [{
                	textStyle:{fontFamily: 'Montserrat'},
                name:"<%=location%>1",
                nameLocation:"middle",//middle start end
                nameGap:20,
                type: "value",
               // data: xData,当type为value时,data不起作用，此时使用的时series中的第一维数据
            }],
            
            
               yAxis: [{
            	   textStyle:{fontFamily: 'Montserrat'},
                type: "value",
                name:"<%=location%>2",
                nameLocation:"middle",//middle start end
                nameGap:20,
            }],
            
            
                series: [{
                //    name: 'price-area',
                    type: 'scatter',
                    symbolSize: 2,
                     itemStyle: {
                         normal: {
                             borderWidth: 0.5,
                             borderColor: '#fff',
                         }
                     },
                    data: location_1
                }]
            };

        // 使用刚指定的配置项和数据显示图表。
        myChart.setOption(option);

        <%//loading%>
        myChart.hideLoading();
</script>

<script type="text/javascript">
        // 基于准备好的dom，初始化echarts实例
        var myChart = echarts.init(document.getElementById('container_right_right'));
        <%//loading%>
        myChart.showLoading({
            text: 'Loading...',
            textStyle: { fontSize : 30 , color: '#444' },
            effectOption: {backgroundColor: 'rgba(0, 0, 0, 0)'}
        })
        option = {
                title: {
                    text: '<%=searchname%>',
                    textStyle:{fontFamily: 'Montserrat'},
                },
                
                 // 工具箱
              toolbox:{
                show:true,
                feature:{
                    // 数据视图
                    dataView:{
                        show:true
                    },
                    // 还原
                    restore:{
                        show:true
                    },
                    // 区域缩放
                    dataZoom:{
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
              
              
                visualMap: {
                	dimension:3,
                    min: exp_min,           //----------------------------
                    max: exp_max,         //----------------------------
                    orient: 'vertical',
                    right: 0,
                    textStyle:{fontFamily: 'Montserrat'},
                    top: 'center',
                    text: ['expression', ''],
                    calculable: true,
                    inRange: {
                        color: ['#ffffff','#002699']
                    },
                },
                
                tooltip: {
                	textStyle:{fontFamily: 'Montserrat'},
                    trigger: 'item',
                    axisPointer: {
                        type: 'cross'
                    },
                    
                 //   formatter:'{c}'
                    
                    formatter: function (params) {
                        //params.value data= [[1,1,'a',5],[1,1,'a',5]]
                        //return 'Cell:'+params.value[2]+'<br> Cluster: '+params.value[3]
                        //data=[{name:'aaa',value:[1,1,3]},{name:'aaa',value:[2,2,2]},],
                        //console.log(params.data.name)
                        //console.log(params.data.value[0])//value:[1,1,3]中第一个数
                        return params.marker+'Cell:'+params.name+'<br>'+params.marker+"Gene expression: "+params.value[3]
                        
                    }
                },
                
                grid: {
                	top:'13%',
                    left: '5%',
                    right: '18%',
                    bottom: '5%',
                    containLabel: true
                  },
                xAxis: [{
                	textStyle:{fontFamily: 'Montserrat'},
                name:"<%=location%>1",
                nameLocation:"middle",//middle start end
                nameGap:20,
                type: "value",
               // data: xData,当type为value时,data不起作用，此时使用的时series中的第一维数据
            }],
            
            
               yAxis: [{
            	   textStyle:{fontFamily: 'Montserrat'},
                type: "value",
                name:"<%=location%>2",
                nameLocation:"middle",//middle start end
                nameGap:20,
            }],
            
            
                series: [{
                //    name: 'price-area',
                    type: 'scatter',
                    symbolSize: 2,
                     itemStyle: {
                         normal: {
                             borderWidth: 0.5,
                             borderColor: '#f2f2f2',
                         }
                     },
                    data: location_1
                }]
            };

        // 使用刚指定的配置项和数据显示图表。
        myChart.setOption(option);

        <%//loading%>
        myChart.hideLoading();
</script>

<script>
	//监听加载状态改变
	document.onreadystatechange = completeLoading;

	//加载状态为complete时移除loading效果
	function completeLoading() {
	    if (document.readyState == "complete") {
	        var loadingMask = document.getElementById('loader-wrapper');
	        loadingMask.parentNode.removeChild(loadingMask);
	    }
	}
	
</script>



</body>
</html>