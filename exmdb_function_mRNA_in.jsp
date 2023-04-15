<%@ page language="java" import="java.io.*,java.sql.*,java.util.*,java.text.DecimalFormat,wp.base.*" contentType="text/html; charset=UTF-8" pageEncoding="UTF-8"%>
<!DOCTYPE html>
<html>
<head>
<meta charset="ISO-8859-1">



<title>Functions</title>

<!-- echarts -->
<script type="text/javascript" language="javascript" src="echart/echarts.min.js"></script>

<style>
body{
	font-family:"Montserrat"
}
</style>

</head>



<body>
<%
//get parameter

String GSEID = "GSE110271";
int i = 0;
String path_GSE = request.getRealPath("/").replace("\\","/");   //get root real path 

if(request.getParameter("GSEID")!=null&&!request.getParameter("GSEID").equals("")){
	GSEID = request.getParameter("GSEID");
}


String Lc = "#cca300";
String Uc = "#002db3";

if(request.getParameter("bottom_color")!=null&&!request.getParameter("bottom_color").equals("")){	
	Uc = request.getParameter("bottom_color");
}
if(request.getParameter("top_color")!=null&&!request.getParameter("top_color").equals("")){	
	Lc = request.getParameter("top_color");
}

// System.out.println("Lc:"+Lc+" Uc:"+Uc);

String searchname="none";
String top = "10";   //默认 显示前10 个功能 
String Sql_volcano = "select * from limma_res where i = \"" + GSEID + "\"";
// System.out.println(Sql_volcano);


if(request.getParameter("top")!=null&&!request.getParameter("top").equals("")){
	top = request.getParameter("top");
	//System.out.println("top:"+top);
}

searchname = "TP53,CDK6,EGFR,VEGFA";

if(request.getParameter("searchname")!=null&&!request.getParameter("searchname").equals("")){	
	searchname = request.getParameter("searchname");
}


StringBuffer sb = new StringBuffer();


int index_update = 0;
if(!searchname.equals("")){
	for (String retval: searchname.split(",")){
		sb.append(dbhello.exmdb_get_entrzID(retval));
		sb.append(",");
	}	
}

String geneids = sb.toString();

String gene_label = "";
String entrid = "";
String gene_sql = "";
String gene_sql_res = "";
gene_label = dbhello.exmdb_get_gene_label(searchname);
if (!gene_label.equals("non_mRNA")){
	entrid = dbhello.exmdb_get_entrzID(searchname);
	gene_sql = "select node2 from ppinet where node1=\""+entrid+"\" union select node1 from ppinet where node2=\""+entrid+"\"";
	geneids = dbString.getSelectInString1(gene_sql);
	
// 	System.out.println("gene_label:"+gene_label);
// 	System.out.println("Entrid:"+entrid);
// 	System.out.println("gene_sql:"+gene_sql);
// 	System.out.println("geneids:"+geneids);
}



// System.out.println(sql);


if(geneids.length()==0){
	geneids = "test";
	out.println("<h2 align=\"center\";}> Sorry, no significant functions enriched!</h2>");
	out.println("<h4 align=\"center\" style='font-weight: 100;color: grey'> The reason for this is that the uploaded RNA is not significantly enriched by the pathway. </h4>");
}

// geneids = geneids.replaceAll("'", "");
// geneids = sb.toString();

String path = request.getRealPath("/").replace("\\","/");   //get root real path 

String annofile_pathway = path+"functionFile/pathwayAnno.txt";
String annofile_go = path+"functionFile/goAnno.txt";
String annofile_hallmark = path+"functionFile/hallmarkAnno.txt";

String annofile_Biocarta = path+"functionFile/Biocarta.txt";
String annofile_BP = path+"functionFile/GO_BP.txt";
String annofile_MF = path+"functionFile/GO_MF.txt";
String annofile_CC = path+"functionFile/GO_CC.txt";
String annofile_Imm = path+"functionFile/Immunesigdb.txt";

//kegg
String values_kegg = "";
String names_kegg = "";

// System.out.println("geneids:" + geneids);
// System.out.println("annofile_pathway: "+ annofile_pathway);


String tableStr_pathway = dbString.getFunctionStr_echarts(geneids, annofile_pathway , Integer.parseInt(top));
if(!tableStr_pathway.equals("@@@")){
values_kegg = tableStr_pathway.split("@@@")[0];
names_kegg = tableStr_pathway.split("@@@")[1];
}


//go
String values_go = "";
String names_go = "";
String tableStr_go = dbString.getFunctionStr_echarts(geneids, annofile_go , Integer.parseInt(top));
if(!tableStr_go.equals("@@@")){
values_go = tableStr_go.split("@@@")[0];
names_go = tableStr_go.split("@@@")[1];
}
	



String values_bp = "";
String names_bp = "";
String tableStr_bp = dbString.getFunctionStr_echarts(geneids, annofile_BP , Integer.parseInt(top));
if(!tableStr_bp.equals("@@@")){
	values_bp = tableStr_bp.split("@@@")[0];
	names_bp = tableStr_bp.split("@@@")[1];
}

String values_mf = "";
String names_mf = "";
String tableStr_mf = dbString.getFunctionStr_echarts(geneids, annofile_MF , Integer.parseInt(top));
if(!tableStr_mf.equals("@@@")){
	values_mf = tableStr_mf.split("@@@")[0];
	names_mf = tableStr_mf.split("@@@")[1];
}

String values_cc = "";
String names_cc = "";
String tableStr_cc = dbString.getFunctionStr_echarts(geneids, annofile_CC , Integer.parseInt(top));
if(!tableStr_cc.equals("@@@")){
	values_cc = tableStr_cc.split("@@@")[0];
	names_cc = tableStr_cc.split("@@@")[1];
}


String values_bio = "";
String names_bio = "";
String tableStr_bio = dbString.getFunctionStr_echarts(geneids, annofile_Biocarta , Integer.parseInt(top));
if(!tableStr_bio.equals("@@@")){
	values_bio = tableStr_bio.split("@@@")[0];
	names_bio = tableStr_bio.split("@@@")[1];
}

String values_Imm = "";
String names_Imm = "";
String tableStr_Imm = dbString.getFunctionStr_echarts(geneids, annofile_Imm , Integer.parseInt(top));
if(!tableStr_bio.equals("@@@")){
	values_Imm = tableStr_Imm.split("@@@")[0];
	names_Imm = tableStr_Imm.split("@@@")[1];
}
	// hallmark
	//ArrayList<String> plist = dbString.getHallmarkPvalues(geneids, annofile_hallmark , 10);
	
	
//}else{
//	out.println(dbString.functionNOfoundStr);
//	out.println("No significantly enriched functions!");
//}



%>


<table style="width:100%">

<tr style="width:100%">
<td style="width:50%">
<!-- go -->

<div id="container_bp" style="min-width:400px;height:400px"></div>
</td>
<td style="width:50%">
<!-- kegg -->
<div id="container_kegg" style="min-width:400px;height:400px"></div>
</td>

</tr>

<tr>
<td style="width:50%">
<div id="container_mf" style="min-width:400px;height:400px"></div>
</td>

<td style="width:50%">
<div id="container_bio" style="min-width:400px;height:400px"></div>
</td>

</tr>

<tr>
<td style="width:50%">
<div id="container_cc" style="min-width:400px;height:400px"></div>
</td>

<td style="width:50%">
<div id="container_Imm" style="min-width:400px;height:400px"></div>
</td>
</tr>

</table>


<!-- go -->      
 <script type="text/javascript">
        // 基于准备好的dom，初始化echarts实例
        var myChart = echarts.init(document.getElementById('container_go'));
option = {
    title: {
        text: ' GeneOntology',
        subtext: '(Top <%=top%> enriched)',
        textStyle:{fontFamily: 'Montserrat',fontSize:15},
    },
    tooltip: {
        trigger: 'axis',
        textStyle:{fontFamily: 'Montserrat',fontSize:15},
        axisPointer: {
            type: 'shadow'
        }
    },
    toolbox: {
        show : true,
        feature : {
            mark : {show: true},
            dataView : {show: true, readOnly: false},
            magicType : {
                show: true,
                type: ['pie', 'funnel']
            },
            restore : {show: true},
            saveAsImage : {show: true}
        }
    },
    grid: {
        left: '2%',
        right: '5%',
        bottom: '12%',
        containLabel: true
    },
    xAxis: {
        type: 'value',
        boundaryGap: [0, 0.01]
    },
    yAxis: {
        type: 'category',
        axisLabel: {
            rotate:  0,  //<%//这个是倾斜角度%>
            interval :0,//　　<%//这里是考虑到x轴文件过多的时候设置的，如果文字太多，默认是间隔显示，设置为0，%>
            fontFamily: 'Montserrat',
            formatter: function (value) {<%// 文字过长用...表示%>
                return (value.length > 40 ? (value.slice(0,40)+"...") : value )
              }
		},
        data: [<%=names_go%>],
        textStyle: {
            fontSize: 1
        }
    },
    visualMap: {
        orient: 'horizontal',
        left: 'center',
        min: 0,
        max: 10,
        textStyle:{fontFamily: 'Montserrat',fontSize:15},
        text: [' ', 'P-Value(-Log10)'],
        // Map the score column to color
        dimension: 0,
        inRange: {
        	color: ['<%=Uc%>', '<%=Lc%>']
        }
      },
    series: [
        {
            name: 'P-value(-log10)',
            type: 'bar',
            barWidth:10,  //住状图宽度
            
             itemStyle: {
                    normal: {
                        color: function(params) {//好，这里就是重头戏了，定义一个list，然后根据所以取得不同的值，这样就实现了，
                            // build a color map as your need.
                            var colorList = [
                                  '#00b359'           
                             // '#C1232B','#B5C334','#FCCE10','#E87C25','#27727B',
                              // '#FE8463','#9BCA63','#FAD860','#F3A43B','#60C0DD',
                              // '#D7504B'//,'#C6E579','#F4E001','#F0805A','#26C0C0'
                            ];
                            return colorList[params.dataIndex]
                        },
　　　　　　　　　　　　　　//以下为是否显示，显示位置和显示格式的设置了
                        
                    }
                },
        
            data: [<%=values_go%>]
        },
        
    ]
};


        // 使用刚指定的配置项和数据显示图表。
        myChart.setOption(option);
    </script>




<!--  kegg       -->
 <script type="text/javascript">
        // 基于准备好的dom，初始化echarts实例
        var myChart = echarts.init(document.getElementById('container_kegg'));
        option = {
        	    title: {
        	        text: 'Kyoto Encyclopedia of Genes and Genomes',
        	        subtext: '(Top <%=top%> enriched)',
        	        textStyle:{fontFamily: 'Montserrat',fontSize:15},
        	    },
        	    tooltip: {
        	        trigger: 'axis',
        	        textStyle:{fontFamily: 'Montserrat',fontSize:15},
        	        axisPointer: {
        	            type: 'shadow'
        	        }
        	    },
        	    toolbox: {
        	        show : true,
        	        feature : {
        	            mark : {show: true},
        	            dataView : {show: true, readOnly: false},
        	            magicType : {
        	                show: true,
        	                type: ['pie', 'funnel']	
        	            },
        	            restore : {show: true},
        	            saveAsImage : {show: true}
        	        }
        	    },
        	    grid: {
        	        left: '5%',
        	        right: '2%',
        	        bottom: '12%',
        	        containLabel: true
        	    },
        	    xAxis: {
        	        type: 'value',
        	        boundaryGap: [0, 0.01]
        	    },
        	    yAxis: {
        	        type: 'category',
        	        
        	        axisLabel: {
	                    rotate:  0,  //<%//这个是倾斜角度%>
	                    interval :0,//　　<%//这里是考虑到x轴文件过多的时候设置的，如果文字太多，默认是间隔显示，设置为0，%>
	                    fontFamily: 'Montserrat',
	                    formatter: function (value) {<%// 文字过长用...表示%>
	                        return (value.length > 40 ? (value.slice(0,40)+"...") : value )
	                      }
	        	},
        	        data: [<%=names_kegg%>]
        	    },
        	    visualMap: {
        	    	textStyle:{fontFamily: 'Montserrat',fontSize:15},
        	        orient: 'horizontal',
        	        left: 'center',
        	        min: 0,
        	        max: 10,
        	        text: [' ', 'P-Value(-Log10)'],
        	        // Map the score column to color
        	        dimension: 0,
        	        inRange: {
        	        	color: ['<%=Uc%>', '<%=Lc%>']
        	        }
        	      },
        	    series: [
        	        {
        	            name: 'P-value(-log10)',
        	            type: 'bar',
        	            barWidth:10,  //住状图宽度
        	            
        	             itemStyle: {
        	                    normal: {
        	                        color: function(params) {//好，这里就是重头戏了，定义一个list，然后根据所以取得不同的值，这样就实现了，
        	                            // build a color map as your need.
        	                            var colorList = [
        	                             '#de6721'
        	                             // '#C1232B','#B5C334','#FCCE10','#E87C25','#27727B',
        	                              // '#FE8463','#9BCA63','#FAD860','#F3A43B','#60C0DD',
        	                              // '#D7504B'//,'#C6E579','#F4E001','#F0805A','#26C0C0'
        	                            ];
        	                            return colorList[params.dataIndex]
        	                        },
        	　　　　　　　　　　　　//以下为是否显示，显示位置和显示格式的设置了
        	                        
        	                    }
        	                },

        	            data: [<%=values_kegg%>]
        	                
        	        },
        	        
        	    ]
        	};



        // 使用刚指定的配置项和数据显示图表。
        myChart.setOption(option);
    </script>

  
 
<script type="text/javascript">
        // 基于准备好的dom，初始化echarts实例
        var myChart = echarts.init(document.getElementById('container_bp'));
        option = {
        	    title: {
        	        text: 'GeneOntology: Biological Process',
        	        textStyle:{fontFamily: 'Montserrat',fontSize:15},
        	        subtext: '(Top <%=top%> enriched)'
        	    },
        	    tooltip: {
        	        trigger: 'axis',
        	        textStyle:{fontFamily: 'Montserrat',fontSize:15},
        	        axisPointer: {
        	            type: 'shadow'
        	        }
        	    },
        	    toolbox: {
        	        show : true,
        	        feature : {
        	            mark : {show: true},
        	            dataView : {show: true, readOnly: false},
        	            magicType : {
        	                show: true,
        	                type: ['pie', 'funnel']	
        	            },
        	            restore : {show: true},
        	            saveAsImage : {show: true}
        	        }
        	    },
//         	    legend: {
//         	        data: ['P-value(-log10)']
//         	    },
        	    grid: {
        	        left: '5%',
        	        right: '2%',
        	        bottom: '12%',
        	        containLabel: true
        	    },
        	    xAxis: {
        	        type: 'value',
        	        boundaryGap: [0, 0.01]
        	    },
        	    yAxis: {
        	        type: 'category',
        	        
        	        axisLabel: {
	                    rotate:  0,  //<%//这个是倾斜角度%>
	                    interval :0,//　　<%//这里是考虑到x轴文件过多的时候设置的，如果文字太多，默认是间隔显示，设置为0，%>
	                    fontFamily: 'Montserrat',
	                    formatter: function (value) {<%// 文字过长用...表示%>
	                        return (value.length > 40 ? (value.slice(0,40)+"...") : value )
	                      }
	        	},
	        	
        	        data: [<%=names_bp%>]
        	    },
        	    visualMap: {
        	    	textStyle:{fontFamily: 'Montserrat',fontSize:15},
        	        orient: 'horizontal',
        	        left: 'center',
        	        min: 0,
        	        max: 10,
        	        text: [' ', 'P-Value(-Log10)'],
        	        // Map the score column to color
        	        dimension: 0,
        	        inRange: {
        	        	color: ['<%=Uc%>', '<%=Lc%>']
        	        }
        	      },
        	    series: [
        	        {
        	            name: 'P-value(-log10)',
        	            type: 'bar',
        	            barWidth:10,  //住状图宽度
        	            
        	             itemStyle: {
        	                    normal: {
        	                        color: function(params) {//好，这里就是重头戏了，定义一个list，然后根据所以取得不同的值，这样就实现了，
        	                            // build a color map as your need.
        	                            var colorList = [
        	                             '#de6721'
        	                             // '#C1232B','#B5C334','#FCCE10','#E87C25','#27727B',
        	                              // '#FE8463','#9BCA63','#FAD860','#F3A43B','#60C0DD',
        	                              // '#D7504B'//,'#C6E579','#F4E001','#F0805A','#26C0C0'
        	                            ];
        	                            return colorList[params.dataIndex]
        	                        },
        	　　　　　　　　　　　　//以下为是否显示，显示位置和显示格式的设置了
        	                        
        	                    }
        	                },

        	            data: [<%=values_bp%>]
        	                
        	        },
        	        
        	    ]
        	};



        // 使用刚指定的配置项和数据显示图表。
        myChart.setOption(option);
    </script>



<script type="text/javascript">
        // 基于准备好的dom，初始化echarts实例
        var myChart = echarts.init(document.getElementById('container_mf'));
        option = {
        	    title: {
        	        text: 'GeneOntology: Molecular Function',
        	        textStyle:{fontFamily: 'Montserrat',fontSize:15},
        	        subtext: '(Top <%=top%> enriched)'
        	    },
        	    tooltip: {
        	        trigger: 'axis',
        	        textStyle:{fontFamily: 'Montserrat',fontSize:15},
        	        axisPointer: {
        	            type: 'shadow'
        	        }
        	    },
        	    toolbox: {
        	        show : true,
        	        feature : {
        	            mark : {show: true},
        	            dataView : {show: true, readOnly: false},
        	            magicType : {
        	                show: true,
        	                type: ['pie', 'funnel']	
        	            },
        	            restore : {show: true},
        	            saveAsImage : {show: true}
        	        }
        	    },
        	    grid: {
        	        left: '5%',
        	        right: '2%',
        	        bottom: '12%',
        	        containLabel: true
        	    },
        	    xAxis: {
        	        type: 'value',
        	        boundaryGap: [0, 0.01]
        	    },
        	    yAxis: {
        	        type: 'category',
        	        axisLabel: {
	                    rotate:  0,  //<%//这个是倾斜角度%>
	                    interval :0,//　　<%//这里是考虑到x轴文件过多的时候设置的，如果文字太多，默认是间隔显示，设置为0，%>
	                    fontFamily: 'Montserrat',
	                    formatter: function (value) {<%// 文字过长用...表示%>
	                        return (value.length > 40 ? (value.slice(0,40)+"...") : value )
	                      }
	        	},
        	        data: [<%=names_mf%>]
        	    },
        	    visualMap: {
        	    	textStyle:{fontFamily: 'Montserrat',fontSize:15},
        	        orient: 'horizontal',
        	        left: 'center',
        	        min: 0,
        	        max: 10,
        	        text: [' ', 'P-Value(-Log10)'],
        	        // Map the score column to color
        	        dimension: 0,
        	        inRange: {
        	          color: ['<%=Uc%>', '<%=Lc%>']
        	        }
        	      },
        	    series: [
        	        {
        	            name: 'P-value(-log10)',
        	            type: 'bar',
        	            barWidth:10,  //住状图宽度
        	            
        	             itemStyle: {
        	                    normal: {
        	                        color: function(params) {//好，这里就是重头戏了，定义一个list，然后根据所以取得不同的值，这样就实现了，
        	                            // build a color map as your need.
        	                            var colorList = [
        	                             '#de6721'
        	                             // '#C1232B','#B5C334','#FCCE10','#E87C25','#27727B',
        	                              // '#FE8463','#9BCA63','#FAD860','#F3A43B','#60C0DD',
        	                              // '#D7504B'//,'#C6E579','#F4E001','#F0805A','#26C0C0'
        	                            ];
        	                            return colorList[params.dataIndex]
        	                        },
        	　　　　　　　　　　　　//以下为是否显示，显示位置和显示格式的设置了
        	                        
        	                    }
        	                },

        	            data: [<%=values_mf%>]
        	                
        	        },
        	        
        	    ]
        	};



        // 使用刚指定的配置项和数据显示图表。
        myChart.setOption(option);
    </script>


<script type="text/javascript">
        // 基于准备好的dom，初始化echarts实例
        var myChart = echarts.init(document.getElementById('container_cc'));
        option = {
        	    title: {
        	        text: 'GeneOntology: Cellular Component',
        	        textStyle:{fontFamily: 'Montserrat',fontSize:15},
        	        subtext: '(Top <%=top%> enriched)'
        	    },
        	    tooltip: {
        	        trigger: 'axis',
        	        textStyle:{fontFamily: 'Montserrat',fontSize:15},
        	        axisPointer: {
        	            type: 'shadow'
        	        }
        	    },
        	    toolbox: {
        	        show : true,
        	        feature : {
        	            mark : {show: true},
        	            dataView : {show: true, readOnly: false},
        	            magicType : {
        	                show: true,
        	                type: ['pie', 'funnel']	
        	            },
        	            restore : {show: true},
        	            saveAsImage : {show: true}
        	        }
        	    },

        	    grid: {
        	        left: '5%',
        	        right: '2%',
        	        bottom: '12%',
        	        containLabel: true
        	    },
        	    xAxis: {
        	        type: 'value',
        	        boundaryGap: [0, 0.01]
        	    },
        	    yAxis: {
        	        type: 'category',
        	        axisLabel: {
	                    rotate:  0,  //<%//这个是倾斜角度%>
	                    interval :0,//　　<%//这里是考虑到x轴文件过多的时候设置的，如果文字太多，默认是间隔显示，设置为0，%>
	                    fontFamily: 'Montserrat',
	                    formatter: function (value) {<%// 文字过长用...表示%>
	                        return (value.length > 40 ? (value.slice(0,40)+"...") : value )
	                      }
	        	},
        	        data: [<%=names_cc%>]
        	    },
        	    visualMap: {
        	    	textStyle:{fontFamily: 'Montserrat',fontSize:15},
        	        orient: 'horizontal',
        	        left: 'center',
        	        min: 0,
        	        max: 10,
        	        text: [' ', 'P-Value(-Log10)'],
        	        // Map the score column to color
        	        dimension: 0,
        	        inRange: {
        	        	color: ['<%=Uc%>', '<%=Lc%>']
        	        }
        	      },
        	    series: [
        	        {
        	            name: 'P-value(-log10)',
        	            type: 'bar',
        	            barWidth:10,  //住状图宽度
        	            
        	             itemStyle: {
        	                    normal: {
        	                    	color: function(params) {
        	                    		return 'rgb(100, 20,'+ ((10-params.value)*20) + ')'
        	                        },
        	　　　　　　　　　　　　//以下为是否显示，显示位置和显示格式的设置了
        	                        
        	                    }
        	                },

        	            data: [<%=values_cc%>]
        	                
        	        },
        	        
        	    ]
        	};



        // 使用刚指定的配置项和数据显示图表。
        myChart.setOption(option);
    </script>




<script type="text/javascript">
        // 基于准备好的dom，初始化echarts实例
        var myChart = echarts.init(document.getElementById('container_bio'));
        option = {
        	    title: {
        	        text: 'Other Pathway [Biocarta,REACTOME,PID,etc.]',
        	        textStyle:{fontFamily: 'Montserrat',fontSize:15},
        	        subtext: '(Top <%=top%> enriched)'
        	    },
        	    tooltip: {
        	        trigger: 'axis',
        	        textStyle:{fontFamily: 'Montserrat',fontSize:15},
        	        axisPointer: {
        	            type: 'shadow'
        	        }
        	    },
        	    toolbox: {
        	        show : true,
        	        feature : {
        	            mark : {show: true},
        	            dataView : {show: true, readOnly: false},
        	            magicType : {
        	                show: true,
        	                type: ['pie', 'funnel']	
        	            },
        	            restore : {show: true},
        	            saveAsImage : {show: true}
        	        }
        	    },

        	    grid: {
        	        left: '5%',
        	        right: '2%',
        	        bottom: '12%',
        	        containLabel: true
        	    },
        	    xAxis: {
        	        type: 'value',
        	        boundaryGap: [0, 0.01]
        	    },
        	    yAxis: {
        	        type: 'category',
        	        axisLabel: {
	                    rotate:  0,  //<%//这个是倾斜角度%>
	                    interval :0,//　　<%//这里是考虑到x轴文件过多的时候设置的，如果文字太多，默认是间隔显示，设置为0，%>
	                    fontFamily: 'Montserrat',
	                    formatter: function (value) {<%// 文字过长用...表示%>
	                        return (value.length > 40 ? (value.slice(0,40)+"...") : value )
	                      }
	        	},
        	        data: [<%=names_bio%>]
        	    },
        	    visualMap: {
        	    	textStyle:{fontFamily: 'Montserrat',fontSize:15},
        	        orient: 'horizontal',
        	        left: 'center',
        	        min: 0,
        	        max: 10,
        	        text: [' ', 'P-Value(-Log10)'],
        	        // Map the score column to color
        	        dimension: 0,
        	        inRange: {
        	        	color: ['<%=Uc%>', '<%=Lc%>']
        	        }
        	      },
        	    series: [
        	        {
        	            name: 'P-value(-log10)',
        	            type: 'bar',
        	            barWidth:10,  //住状图宽度
        	            
        	             itemStyle: {
        	                    normal: {
        	                        color: function(params) {//好，这里就是重头戏了，定义一个list，然后根据所以取得不同的值，这样就实现了，
        	                            // build a color map as your need.
        	                            var colorList = [
        	                             '#de6721'
        	                             // '#C1232B','#B5C334','#FCCE10','#E87C25','#27727B',
        	                              // '#FE8463','#9BCA63','#FAD860','#F3A43B','#60C0DD',
        	                              // '#D7504B'//,'#C6E579','#F4E001','#F0805A','#26C0C0'
        	                            ];
        	                            return colorList[params.dataIndex]
        	                        },
        	　　　　　　　　　　　　//以下为是否显示，显示位置和显示格式的设置了
        	                        
        	                    }
        	                },

        	            data: [<%=values_bio%>]
        	                
        	        },
        	        
        	    ]
        	};



        // 使用刚指定的配置项和数据显示图表。
        myChart.setOption(option);
        myChart.on('click', function (params) {
			console.log(params.data);
        })
    </script>

<script type="text/javascript">
        // 基于准备好的dom，初始化echarts实例
        var myChart = echarts.init(document.getElementById('container_Imm'));
        option = {
        	    title: {
        	        text: 'Msigdb : Immune pathway',
        	        textStyle:{fontFamily: 'Montserrat',fontSize:15},
        	        subtext: '(Top <%=top%> enriched)'
        	    },
        	    tooltip: {
        	        trigger: 'axis',
        	        textStyle:{fontFamily: 'Montserrat',fontSize:15},
        	        axisPointer: {
        	            type: 'shadow'
        	        }
        	    },
        	    toolbox: {
        	        show : true,
        	        feature : {
        	            mark : {show: true},
        	            dataView : {show: true, readOnly: false},
        	            magicType : {
        	                show: true,
        	                type: ['pie', 'funnel']	
        	            },
        	            restore : {show: true},
        	            saveAsImage : {show: true}
        	        }
        	    },
        	    grid: {
        	        left: '5%',
        	        right: '2%',
        	        bottom: '12%',
        	        containLabel: true
        	    },
        	    xAxis: {
        	        type: 'value',
        	        boundaryGap: [0, 0.01]
        	    },
        	    yAxis: {
        	        type: 'category',
        	        axisLabel: {
	                    rotate:  0,  //<%//这个是倾斜角度%>
	                    interval :0,//　　<%//这里是考虑到x轴文件过多的时候设置的，如果文字太多，默认是间隔显示，设置为0，%>
	                    fontFamily: 'Montserrat',
	                    formatter: function (value) {<%// 文字过长用...表示%>
	                        return (value.length > 40 ? (value.slice(0,40)+"...") : value )
	                      }
	        	},
        	        data: [<%=names_Imm%>]
        	    },
        	    visualMap: {
        	    	textStyle:{fontFamily: 'Montserrat',fontSize:15},
        	        orient: 'horizontal',
        	        left: 'center',
        	        min: 0,
        	        max: 10,
        	        text: [' ', 'P-Value(-Log10)'],
        	        // Map the score column to color
        	        dimension: 0,
        	        inRange: {
        	          color: ['<%=Uc%>', '<%=Lc%>']
        	        }
        	      },
        	    series: [
        	        {
        	            name: 'P-value(-log10)',
        	            type: 'bar',
        	            barWidth:10,  //住状图宽度
        	            
        	             itemStyle: {
        	                    normal: {
        	                        color: function(params) {//好，这里就是重头戏了，定义一个list，然后根据所以取得不同的值，这样就实现了，
        	                            // build a color map as your need.
        	                            var colorList = [
        	                             '#de6721'
        	                             // '#C1232B','#B5C334','#FCCE10','#E87C25','#27727B',
        	                              // '#FE8463','#9BCA63','#FAD860','#F3A43B','#60C0DD',
        	                              // '#D7504B'//,'#C6E579','#F4E001','#F0805A','#26C0C0'
        	                            ];
        	                            return colorList[params.dataIndex]
        	                        },
        	　　　　　　　　　　　　//以下为是否显示，显示位置和显示格式的设置了
        	                        
        	                    }
        	                },

        	            data: [<%=values_Imm%>]
        	                
        	        },
        	        
        	    ]
        	};



        // 使用刚指定的配置项和数据显示图表。
        myChart.setOption(option);
    </script>
</body>
</html>