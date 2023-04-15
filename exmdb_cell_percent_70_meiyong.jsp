<%@ page language="java" import="java.io.*,java.sql.*,java.util.*,java.text.DecimalFormat,wp.base.*" contentType="text/html; charset=UTF-8" pageEncoding="UTF-8"%>
<!doctype html>
<html lang="zxx" style="height: 100%">
    <head>
        <meta charset="utf-8">
    </head>
    <body style="height: 100%; margin: 0">
        <div id="container" style="height: 100%"></div>
        <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/echarts@5/dist/echarts.min.js"></script>
        <!-- Uncomment this line if you want to dataTool extension
        <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/echarts@5/dist/extension/dataTool.min.js"></script>
        -->
        <!-- Uncomment this line if you want to use gl extension
        <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/echarts-gl@2/dist/echarts-gl.min.js"></script>
        -->
        <!-- Uncomment this line if you want to echarts-stat extension
        <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/echarts-stat@latest/dist/ecStat.min.js"></script>
        -->
        <!-- Uncomment this line if you want to use map
        <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/echarts@5/map/js/china.js"></script>
        <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/echarts@5/map/js/world.js"></script>
        -->
        <!-- Uncomment these two lines if you want to use bmap extension
        <script type="text/javascript" src="https://api.map.baidu.com/api?v=2.0&ak=<Your Key Here>"></script>
        <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/echarts@{{version}}/dist/extension/bmap.min.js"></script>
        -->
<% 


String qht = dbhello.getpercent("select * from mean_exp_cell where genesymbol = \"PLIN2: AEL_GSE142213\"");


%>
        <script type="text/javascript">
			var dom = document.getElementById("container");
			var myChart = echarts.init(dom);
			var app = {};
			
			var option;
			
			
			
			option = {
			title: {text:"Gene: QIHAITAO",textStyle:{fontFamily: 'Montserrat',fontSize:20}},
			legend: {
					textStyle:{fontFamily: 'Montserrat',fontSize:15},
					left: 'left',
					bottom: 0,
			},
			  tooltip: {textStyle:{fontFamily: 'Montserrat',fontSize:15}},
			  dataset: {
			    source: [
			      ['Cell Type', 'Malignant_cell','Mono_Macro_cell','NK_cell','Plasma_cell','B_cell','CD4Tconv_cell','CD8T_cell','Erythrocytes_cell','Tprolif_cell','EryPro_cell','GMP_cell','HSC_cell','Progenitor_cell','Promonocyte_cell','CD8Tex_cell','DC_cell','Endothelial_cell','Fibroblasts_cell','Mast_cell','Melanocytes_cell','Myofibroblasts_cell','Treg_cell','Epithelial_cell','Neutrophils_cell','Oligodendrocyte_cell','Hepatic_progenitor_cell','AC_like_Malignant_cell','OC_like_Malignant_cell','OPC_like_Malignant_cell','NB_like_Malignant_cell','Neuron_cell','MES_like_Malignant_cell','NPC_like_Malignant_cell','Microglia_cell','Astrocyte_cell','OPC_cell','Vascular_cell','Myocyte_cell','ILC_cell','Others_cell','Keratinocytes_cell','Pericytes_cell','Secretory_glandular_cell','SMC_cell','Alveolar_cell','Acinar_cell','Ductal_cell','Gland_mucous_cell','Pit_mucous_cell','Stromal_cell','Myeloids_cell','T_cell','iCAFs_cell','myCAFs_cell','Epithelial_Basal_cell','T_unassigned_cell','NKT_cell','T_Cycling_cell','CD4','Tfh_cell','dPVL_cell','imPVL_cell','Myoepithelial_cell','Epithelial_Basal_Cycling_cell','Epithelial_Luminal_Mature_cell','Mesenchymal_stem_cell','Megakaryocytes','Luminal_progenitor_cell','T_helper_cell','Luminal_cell'],
					<%=qht%>
			    ]
			  },
			  color: ["#FDF5B0","#3D6DFB"],
			  series: [
				    {type: 'pie',radius: '12%',center: ['5%', '10%'],encode: {itemName: 'Cell Type',value:'Malignant_cell'},name:'Malignant_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['15%', '10%'],encode: {itemName: 'Cell Type',value:'Mono_Macro_cell'},name:'Mono_Macro_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['25%', '10%'],encode: {itemName: 'Cell Type',value:'NK_cell'},name:'NK_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['35%', '10%'],encode: {itemName: 'Cell Type',value:'Plasma_cell'},name:'Plasma_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['45%', '10%'],encode: {itemName: 'Cell Type',value:'B_cell'},name:'B_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['55%', '10%'],encode: {itemName: 'Cell Type',value:'CD4Tconv_cell'},name:'CD4Tconv_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['65%', '10%'],encode: {itemName: 'Cell Type',value:'CD8T_cell'},name:'CD8T_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['75%', '10%'],encode: {itemName: 'Cell Type',value:'Erythrocytes_cell'},name:'Erythrocytes_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['85%', '10%'],encode: {itemName: 'Cell Type',value:'Tprolif_cell'},name:'Tprolif_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['95%', '10%'],encode: {itemName: 'Cell Type',value:'EryPro_cell'},name:'EryPro_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['5%', '23%'],encode: {itemName: 'Cell Type',value:'GMP_cell'},name:'GMP_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['15%', '23%'],encode: {itemName: 'Cell Type',value:'HSC_cell'},name:'HSC_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['25%', '23%'],encode: {itemName: 'Cell Type',value:'Progenitor_cell'},name:'Progenitor_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['35%', '23%'],encode: {itemName: 'Cell Type',value:'Promonocyte_cell'},name:'Promonocyte_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['45%', '23%'],encode: {itemName: 'Cell Type',value:'CD8Tex_cell'},name:'CD8Tex_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['55%', '23%'],encode: {itemName: 'Cell Type',value:'DC_cell'},name:'DC_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['65%', '23%'],encode: {itemName: 'Cell Type',value:'Endothelial_cell'},name:'Endothelial_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['75%', '23%'],encode: {itemName: 'Cell Type',value:'Fibroblasts_cell'},name:'Fibroblasts_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['85%', '23%'],encode: {itemName: 'Cell Type',value:'Mast_cell'},name:'Mast_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['95%', '23%'],encode: {itemName: 'Cell Type',value:'Melanocytes_cell'},name:'Melanocytes_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['5%', '36%'],encode: {itemName: 'Cell Type',value:'Myofibroblasts_cell'},name:'Myofibroblasts_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['15%', '36%'],encode: {itemName: 'Cell Type',value:'Treg_cell'},name:'Treg_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['25%', '36%'],encode: {itemName: 'Cell Type',value:'Epithelial_cell'},name:'Epithelial_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['35%', '36%'],encode: {itemName: 'Cell Type',value:'Neutrophils_cell'},name:'Neutrophils_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['45%', '36%'],encode: {itemName: 'Cell Type',value:'Oligodendrocyte_cell'},name:'Oligodendrocyte_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['55%', '36%'],encode: {itemName: 'Cell Type',value:'Hepatic_progenitor_cell'},name:'Hepatic_progenitor_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['65%', '36%'],encode: {itemName: 'Cell Type',value:'AC_like_Malignant_cell'},name:'AC_like_Malignant_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['75%', '36%'],encode: {itemName: 'Cell Type',value:'OC_like_Malignant_cell'},name:'OC_like_Malignant_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['85%', '36%'],encode: {itemName: 'Cell Type',value:'OPC_like_Malignant_cell'},name:'OPC_like_Malignant_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['95%', '36%'],encode: {itemName: 'Cell Type',value:'NB_like_Malignant_cell'},name:'NB_like_Malignant_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['5%', '49%'],encode: {itemName: 'Cell Type',value:'Neuron_cell'},name:'Neuron_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['15%', '49%'],encode: {itemName: 'Cell Type',value:'MES_like_Malignant_cell'},name:'MES_like_Malignant_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['25%', '49%'],encode: {itemName: 'Cell Type',value:'NPC_like_Malignant_cell'},name:'NPC_like_Malignant_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['35%', '49%'],encode: {itemName: 'Cell Type',value:'Microglia_cell'},name:'Microglia_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['45%', '49%'],encode: {itemName: 'Cell Type',value:'Astrocyte_cell'},name:'Astrocyte_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['55%', '49%'],encode: {itemName: 'Cell Type',value:'OPC_cell'},name:'OPC_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['65%', '49%'],encode: {itemName: 'Cell Type',value:'Vascular_cell'},name:'Vascular_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['75%', '49%'],encode: {itemName: 'Cell Type',value:'Myocyte_cell'},name:'Myocyte_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['85%', '49%'],encode: {itemName: 'Cell Type',value:'ILC_cell'},name:'ILC_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['95%', '49%'],encode: {itemName: 'Cell Type',value:'Others_cell'},name:'Others_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['5%', '62%'],encode: {itemName: 'Cell Type',value:'Keratinocytes_cell'},name:'Keratinocytes_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['15%', '62%'],encode: {itemName: 'Cell Type',value:'Pericytes_cell'},name:'Pericytes_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['25%', '62%'],encode: {itemName: 'Cell Type',value:'Secretory_glandular_cell'},name:'Secretory_glandular_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['35%', '62%'],encode: {itemName: 'Cell Type',value:'SMC_cell'},name:'SMC_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['45%', '62%'],encode: {itemName: 'Cell Type',value:'Alveolar_cell'},name:'Alveolar_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['55%', '62%'],encode: {itemName: 'Cell Type',value:'Acinar_cell'},name:'Acinar_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['65%', '62%'],encode: {itemName: 'Cell Type',value:'Ductal_cell'},name:'Ductal_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['75%', '62%'],encode: {itemName: 'Cell Type',value:'Gland_mucous_cell'},name:'Gland_mucous_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['85%', '62%'],encode: {itemName: 'Cell Type',value:'Pit_mucous_cell'},name:'Pit_mucous_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['95%', '62%'],encode: {itemName: 'Cell Type',value:'Stromal_cell'},name:'Stromal_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['5%', '75%'],encode: {itemName: 'Cell Type',value:'Myeloids_cell'},name:'Myeloids_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['15%', '75%'],encode: {itemName: 'Cell Type',value:'T_cell'},name:'T_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['25%', '75%'],encode: {itemName: 'Cell Type',value:'iCAFs_cell'},name:'iCAFs_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['35%', '75%'],encode: {itemName: 'Cell Type',value:'myCAFs_cell'},name:'myCAFs_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['45%', '75%'],encode: {itemName: 'Cell Type',value:'Epithelial_Basal_cell'},name:'Epithelial_Basal_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['55%', '75%'],encode: {itemName: 'Cell Type',value:'T_unassigned_cell'},name:'T_unassigned_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['65%', '75%'],encode: {itemName: 'Cell Type',value:'NKT_cell'},name:'NKT_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['75%', '75%'],encode: {itemName: 'Cell Type',value:'T_Cycling_cell'},name:'T_Cycling_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['85%', '75%'],encode: {itemName: 'Cell Type',value:'CD4'},name:'CD4',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['95%', '75%'],encode: {itemName: 'Cell Type',value:'Tfh_cell'},name:'Tfh_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['5%', '88%'],encode: {itemName: 'Cell Type',value:'dPVL_cell'},name:'dPVL_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['15%', '88%'],encode: {itemName: 'Cell Type',value:'imPVL_cell'},name:'imPVL_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['25%', '88%'],encode: {itemName: 'Cell Type',value:'Myoepithelial_cell'},name:'Myoepithelial_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['35%', '88%'],encode: {itemName: 'Cell Type',value:'Epithelial_Basal_Cycling_cell'},name:'Epithelial_Basal_Cycling_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['45%', '88%'],encode: {itemName: 'Cell Type',value:'Epithelial_Luminal_Mature_cell'},name:'Epithelial_Luminal_Mature_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['55%', '88%'],encode: {itemName: 'Cell Type',value:'Mesenchymal_stem_cell'},name:'Mesenchymal_stem_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['65%', '88%'],encode: {itemName: 'Cell Type',value:'Megakaryocytes'},name:'Megakaryocytes',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['75%', '88%'],encode: {itemName: 'Cell Type',value:'Luminal_progenitor_cell'},name:'Luminal_progenitor_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['85%', '88%'],encode: {itemName: 'Cell Type',value:'T_helper_cell'},name:'T_helper_cell',label: { show: false}},
				    {type: 'pie',radius: '12%',center: ['95%', '88%'],encode: {itemName: 'Cell Type',value:'Luminal_cell'},name:'Luminal_cell',label: { show: false}},
			  ]
			};
			
			if (option && typeof option === 'object') {
			    myChart.setOption(option);
			}

        </script>
    </body>
</html>
    