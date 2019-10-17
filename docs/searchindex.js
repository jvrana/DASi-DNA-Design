Search.setIndex({docnames:["_autosummary/dasi.design","_autosummary/dasi.utils","command_line","constants","cost","design","developer/api_reference","developer/changelog","exceptions","generated/dasi.cost.params","generated/dasi.cost.utils","generated/dasi.design.design_algorithms","generated/dasi.design.graph_builder","generated/dasi.design.optimize","generated/dasi.design.plotter","generated/dasi.models.AlignmentContainer","generated/dasi.models.AlignmentContainerFactory","generated/dasi.models.AlignmentGroup","generated/dasi.models.AlignmentGroupBase","generated/dasi.models.Assembly","generated/dasi.models.AssemblyNode","generated/dasi.models.Molecule","generated/dasi.models.MoleculeType","generated/dasi.models.MultiPCRProductAlignmentGroup","generated/dasi.models.PCRProductAlignmentGroup","generated/dasi.models.Reaction","generated/dasi.models.alignment","generated/dasi.models.alignment_container","generated/dasi.utils.networkx","generated/dasi.utils.networkx.exceptions","generated/dasi.utils.networkx.shortest_path","generated/dasi.utils.networkx.utils","generated/dasi.utils.npdf","generated/dasi.utils.region","generated/dasi.utils.sequence_design","guidelines","index","models","schemas/schemas","usage","utils"],envversion:{"sphinx.domains.c":1,"sphinx.domains.changeset":1,"sphinx.domains.citation":1,"sphinx.domains.cpp":1,"sphinx.domains.javascript":1,"sphinx.domains.math":2,"sphinx.domains.python":1,"sphinx.domains.rst":1,"sphinx.domains.std":1,"sphinx.ext.intersphinx":1,"sphinx.ext.viewcode":1,nbsphinx:1,sphinx:56},filenames:["_autosummary/dasi.design.rst","_autosummary/dasi.utils.rst","command_line.rst","constants.rst","cost.rst","design.rst","developer/api_reference.rst","developer/changelog.md","exceptions.rst","generated/dasi.cost.params.rst","generated/dasi.cost.utils.rst","generated/dasi.design.design_algorithms.rst","generated/dasi.design.graph_builder.rst","generated/dasi.design.optimize.rst","generated/dasi.design.plotter.rst","generated/dasi.models.AlignmentContainer.rst","generated/dasi.models.AlignmentContainerFactory.rst","generated/dasi.models.AlignmentGroup.rst","generated/dasi.models.AlignmentGroupBase.rst","generated/dasi.models.Assembly.rst","generated/dasi.models.AssemblyNode.rst","generated/dasi.models.Molecule.rst","generated/dasi.models.MoleculeType.rst","generated/dasi.models.MultiPCRProductAlignmentGroup.rst","generated/dasi.models.PCRProductAlignmentGroup.rst","generated/dasi.models.Reaction.rst","generated/dasi.models.alignment.rst","generated/dasi.models.alignment_container.rst","generated/dasi.utils.networkx.rst","generated/dasi.utils.networkx.exceptions.rst","generated/dasi.utils.networkx.shortest_path.rst","generated/dasi.utils.networkx.utils.rst","generated/dasi.utils.npdf.rst","generated/dasi.utils.region.rst","generated/dasi.utils.sequence_design.rst","guidelines.rst","index.rst","models.rst","schemas/schemas.rst","usage.rst","utils.rst"],objects:{"dasi.command_line":{DasiCLI:[2,1,1,""]},"dasi.command_line.DasiCLI":{_get_span_cost:[2,2,1,""],run:[2,2,1,""],version:[2,2,1,""]},"dasi.constants":{Constants:[3,1,1,""]},"dasi.constants.Constants":{FRAGMENT:[3,3,1,""],GAP:[3,3,1,""],INF:[3,3,1,""],MAX_HOMOLOGY:[3,3,1,""],MIN_OVERLAP:[3,3,1,""],MISSING:[3,3,1,""],OVERLAP:[3,3,1,""],PCR_PRODUCT:[3,3,1,""],PCR_PRODUCT_WITH_LEFT_PRIMER:[3,3,1,""],PCR_PRODUCT_WITH_PRIMERS:[3,3,1,""],PCR_PRODUCT_WITH_RIGHT_PRIMER:[3,3,1,""],PRIMER:[3,3,1,""],PRIMER_EXTENSION_PRODUCT:[3,3,1,""],PRIMER_EXTENSION_PRODUCT_WITH_LEFT_PRIMER:[3,3,1,""],PRIMER_EXTENSION_PRODUCT_WITH_RIGHT_PRIMER:[3,3,1,""],PRIMER_MIN_BIND:[3,3,1,""],SHARED_FRAGMENT:[3,3,1,""],TEMPLATE:[3,3,1,""]},"dasi.cost":{CostBuilder:[4,1,1,""],PrimerCostModel:[4,1,1,""],SpanCost:[4,1,1,""],SynthesisCostModel:[4,1,1,""],decoder:[4,4,1,""],encoder:[4,4,1,""],utils:[10,0,0,"-"]},"dasi.cost.CostBuilder":{compute:[4,2,1,""],cost:[4,2,1,""],from_json:[4,2,1,""],plot:[4,2,1,""],to_df:[4,2,1,""]},"dasi.cost.PrimerCostModel":{compute:[4,2,1,""],cost:[4,2,1,""],from_json:[4,2,1,""],plot:[4,2,1,""],to_df:[4,2,1,""]},"dasi.cost.SpanCost":{compute:[4,2,1,""],cost:[4,2,1,""],from_json:[4,2,1,""],plot:[4,2,1,""],to_df:[4,2,1,""]},"dasi.cost.SynthesisCostModel":{compute:[4,2,1,""],cost:[4,2,1,""],from_json:[4,2,1,""],plot:[4,2,1,""],to_df:[4,2,1,""]},"dasi.cost.utils":{df_to_np_ranged:[10,4,1,""],duplicates:[10,4,1,""],flatten_axis:[10,4,1,""],square_broadcast:[10,4,1,""]},"dasi.design":{Design:[5,1,1,""],DesignResult:[5,1,1,""],LibraryDesign:[5,1,1,""],design_algorithms:[11,0,0,"-"],graph_builder:[12,0,0,"-"],optimize:[13,0,0,"-"],plotter:[14,0,0,"-"]},"dasi.design.Design":{_assemble_graphs_without_threads:[5,2,1,""],_optimize_without_threads:[5,2,1,""],_results:[5,3,1,""],_seqdb:[5,3,1,""],add_fragments:[5,2,1,""],add_primers:[5,2,1,""],add_queries:[5,2,1,""],add_templates:[5,2,1,""],compile:[5,2,1,""],container_list:[5,2,1,""],containers:[5,2,1,""],filter_linear_records:[5,2,1,""],filter_perfect_subject:[5,2,1,""],graphs:[5,3,1,""],n_jobs:[5,3,1,""],query_keys:[5,2,1,""],span_cost:[5,3,1,""]},"dasi.design.DesignResult":{add_assemblies:[5,2,1,""],add_assembly:[5,2,1,""],assemblies:[5,2,1,""]},"dasi.design.LibraryDesign":{_assemble_graphs_without_threads:[5,2,1,""],_get_iter_non_repeats:[5,2,1,""],_optimize_without_threads:[5,2,1,""],_share_query_blast:[5,2,1,""],add_fragments:[5,2,1,""],add_primers:[5,2,1,""],add_queries:[5,2,1,""],add_templates:[5,2,1,""],compile:[5,2,1,""],compile_library:[5,2,1,""],container_list:[5,2,1,""],containers:[5,2,1,""],filter_linear_records:[5,2,1,""],filter_perfect_subject:[5,2,1,""],optimize_library:[5,2,1,""],query_keys:[5,2,1,""]},"dasi.design.design_algorithms":{assemble_graph:[11,4,1,""],multiprocessing_assemble_graph:[11,4,1,""],multiprocessing_optimize_graph:[11,4,1,""]},"dasi.design.graph_builder":{AssemblyGraphBuilder:[12,1,1,""]},"dasi.design.graph_builder.AssemblyGraphBuilder":{_batch_add_edge_costs:[12,2,1,""],add_edge:[12,2,1,""],add_node:[12,2,1,""]},"dasi.design.optimize":{_check_paths:[13,4,1,""],_multinode_to_shortest_path:[13,4,1,""],_nodes_to_fullpaths:[13,4,1,""],index_slice:[13,4,1,""]},"dasi.exceptions":{AlignmentContainerException:[8,5,1,""],AlignmentException:[8,5,1,""],DASiException:[8,5,1,""],DASiWarning:[8,5,1,""],DasiCostParameterValidationError:[8,5,1,""],DasiDesignException:[8,5,1,""],DasiInvalidMolecularAssembly:[8,5,1,""],DasiNoPrimerPairsException:[8,5,1,""],DasiSequenceDesignException:[8,5,1,""]},"dasi.exceptions.AlignmentContainerException":{with_traceback:[8,2,1,""]},"dasi.exceptions.AlignmentException":{with_traceback:[8,2,1,""]},"dasi.exceptions.DASiException":{with_traceback:[8,2,1,""]},"dasi.exceptions.DASiWarning":{with_traceback:[8,2,1,""]},"dasi.exceptions.DasiCostParameterValidationError":{with_traceback:[8,2,1,""]},"dasi.exceptions.DasiDesignException":{with_traceback:[8,2,1,""]},"dasi.exceptions.DasiInvalidMolecularAssembly":{with_traceback:[8,2,1,""]},"dasi.exceptions.DasiNoPrimerPairsException":{with_traceback:[8,2,1,""]},"dasi.exceptions.DasiSequenceDesignException":{with_traceback:[8,2,1,""]},"dasi.models":{AlignmentContainer:[15,1,1,""],AlignmentContainerFactory:[16,1,1,""],AlignmentGroup:[17,1,1,""],AlignmentGroupBase:[18,1,1,""],Assembly:[19,1,1,""],AssemblyNode:[20,1,1,""],Molecule:[21,1,1,""],MoleculeType:[22,1,1,""],MultiPCRProductAlignmentGroup:[23,1,1,""],PCRProductAlignmentGroup:[24,1,1,""],Reaction:[25,1,1,""],alignment:[26,0,0,"-"],alignment_container:[27,0,0,"-"]},"dasi.models.AlignmentContainer":{__init__:[15,2,1,""],_alignment_hash:[15,2,1,""],_create_pcr_product_alignment:[15,2,1,""],expand:[15,2,1,""],expand_overlaps:[15,2,1,""],expand_primer_pairs:[15,2,1,""],find_groups_by_pos:[15,2,1,""],freeze:[15,2,1,""],get_groups_by_types:[15,2,1,""],groups_by_type:[15,2,1,""],redundent_alignment_groups:[15,2,1,""],types:[15,2,1,""],unfreeze:[15,2,1,""]},"dasi.models.AlignmentContainerFactory":{__init__:[16,2,1,""],alignments:[16,2,1,""],containers:[16,2,1,""],load_blast_json:[16,2,1,""]},"dasi.models.AlignmentGroup":{__init__:[17,2,1,""],query_key:[17,2,1,""],sub_region:[17,2,1,""],subject_keys:[17,2,1,""],subject_regions:[17,2,1,""]},"dasi.models.AlignmentGroupBase":{__init__:[18,2,1,""],query_key:[18,2,1,""],sub_region:[18,2,1,""],subject_keys:[18,2,1,""],subject_regions:[18,2,1,""]},"dasi.models.Assembly":{__init__:[19,2,1,""],_head:[19,2,1,""]},"dasi.models.AssemblyNode":{__init__:[20,2,1,""],_asdict:[20,2,1,""],_make:[20,2,1,""],_replace:[20,2,1,""],count:[20,2,1,""],expandable:[20,2,1,""],index:[20,2,1,""],overhang:[20,2,1,""],type:[20,2,1,""]},"dasi.models.Molecule":{__init__:[21,2,1,""]},"dasi.models.MoleculeType":{__init__:[22,2,1,""]},"dasi.models.MultiPCRProductAlignmentGroup":{__init__:[23,2,1,""],query_key:[23,2,1,""],sub_region:[23,2,1,""],subject_keys:[23,2,1,""],subject_regions:[23,2,1,""]},"dasi.models.PCRProductAlignmentGroup":{__init__:[24,2,1,""],query_key:[24,2,1,""],sub_region:[24,2,1,""],subject_keys:[24,2,1,""],subject_regions:[24,2,1,""]},"dasi.models.Reaction":{__init__:[25,2,1,""],inputs:[25,3,1,""],outputs:[25,3,1,""]},"dasi.models.alignment":{Alignment:[26,1,1,""],AlignmentGroup:[26,1,1,""],AlignmentGroupBase:[26,1,1,""],MultiPCRProductAlignmentGroup:[26,1,1,""],PCRProductAlignmentGroup:[26,1,1,""],RepresentsMolecule:[26,1,1,""]},"dasi.models.alignment.Alignment":{copy:[26,2,1,""],sub_region:[26,2,1,""]},"dasi.models.alignment.AlignmentGroup":{query_key:[26,2,1,""],sub_region:[26,2,1,""],subject_keys:[26,2,1,""],subject_regions:[26,2,1,""]},"dasi.models.alignment.AlignmentGroupBase":{query_key:[26,2,1,""],sub_region:[26,2,1,""],subject_keys:[26,2,1,""],subject_regions:[26,2,1,""]},"dasi.models.alignment.MultiPCRProductAlignmentGroup":{query_key:[26,2,1,""],sub_region:[26,2,1,""],subject_keys:[26,2,1,""],subject_regions:[26,2,1,""]},"dasi.models.alignment.PCRProductAlignmentGroup":{query_key:[26,2,1,""],sub_region:[26,2,1,""],subject_keys:[26,2,1,""],subject_regions:[26,2,1,""]},"dasi.models.alignment_container":{AlignmentContainer:[27,1,1,""],AlignmentContainerFactory:[27,1,1,""],blast_to_region:[27,4,1,""]},"dasi.models.alignment_container.AlignmentContainer":{_alignment_hash:[27,2,1,""],_create_pcr_product_alignment:[27,2,1,""],expand:[27,2,1,""],expand_overlaps:[27,2,1,""],expand_primer_pairs:[27,2,1,""],find_groups_by_pos:[27,2,1,""],freeze:[27,2,1,""],get_groups_by_types:[27,2,1,""],groups_by_type:[27,2,1,""],redundent_alignment_groups:[27,2,1,""],types:[27,2,1,""],unfreeze:[27,2,1,""]},"dasi.models.alignment_container.AlignmentContainerFactory":{alignments:[27,2,1,""],containers:[27,2,1,""],load_blast_json:[27,2,1,""]},"dasi.utils":{bisect_between:[40,4,1,""],bisect_slice_between:[40,4,1,""],networkx:[28,0,0,"-"],npdf:[32,0,0,"-"],perfect_subject:[40,4,1,""],region:[33,0,0,"-"],sequence_design:[34,0,0,"-"],sort_with_keys:[40,4,1,""]},"dasi.utils.networkx":{exceptions:[29,0,0,"-"],shortest_path:[30,0,0,"-"],utils:[31,0,0,"-"]},"dasi.utils.networkx.exceptions":{NetworkxUtilsException:[29,5,1,""]},"dasi.utils.networkx.exceptions.NetworkxUtilsException":{with_traceback:[29,2,1,""]},"dasi.utils.networkx.shortest_path":{multipoint_shortest_path:[30,4,1,""],sympy_dijkstras:[30,4,1,""]},"dasi.utils.networkx.utils":{divide:[31,4,1,""],find_all_min_paths:[31,4,1,""],select_from_arrs:[31,4,1,""],sort_cycle:[31,4,1,""]},"dasi.utils.npdf":{Null:[32,1,1,""],NumpyDataFrame:[32,1,1,""],NumpyDataFrameException:[32,5,1,""],NumpyDataFrameIndexer:[32,1,1,""]},"dasi.utils.npdf.NumpyDataFrame":{aggregate:[32,2,1,""],append:[32,2,1,""],apply:[32,2,1,""],apply_to_col_names:[32,2,1,""],col:[32,2,1,""],columns:[32,2,1,""],concat:[32,2,1,""],copy:[32,2,1,""],data:[32,2,1,""],dump:[32,2,1,""],dumps:[32,2,1,""],fill_value:[32,2,1,""],get:[32,2,1,""],group_apply:[32,2,1,""],hstack:[32,2,1,""],items:[32,2,1,""],keys:[32,2,1,""],load:[32,2,1,""],loads:[32,2,1,""],merge:[32,2,1,""],prefix:[32,2,1,""],reshape:[32,2,1,""],shape:[32,2,1,""],stack:[32,2,1,""],suffix:[32,2,1,""],to_df:[32,2,1,""],update:[32,2,1,""],validate:[32,2,1,""],values:[32,2,1,""],vstack:[32,2,1,""]},"dasi.utils.npdf.NumpyDataFrameException":{with_traceback:[32,2,1,""]},"dasi.utils.npdf.NumpyDataFrameIndexer":{get:[32,2,1,""],items:[32,2,1,""],keys:[32,2,1,""],values:[32,2,1,""]},"dasi.utils.region":{EmptySpan:[33,1,1,""],Region:[33,1,1,""],Span:[33,1,1,""],SpanError:[33,5,1,""]},"dasi.utils.region.EmptySpan":{"new":[33,2,1,""],a:[33,2,1,""],b:[33,2,1,""],bounds:[33,2,1,""],c:[33,2,1,""],connecting_span:[33,2,1,""],consecutive:[33,2,1,""],context_length:[33,2,1,""],cyclic:[33,2,1,""],differences:[33,2,1,""],force_context:[33,2,1,""],get_slice:[33,2,1,""],get_slice_iter:[33,2,1,""],i:[33,2,1,""],index:[33,2,1,""],intersection:[33,2,1,""],invert:[33,2,1,""],overlaps_with:[33,2,1,""],ranges:[33,2,1,""],reindex:[33,2,1,""],same_context:[33,2,1,""],slices:[33,2,1,""],sub:[33,2,1,""],t:[33,2,1,""]},"dasi.utils.region.Region":{"new":[33,2,1,""],a:[33,2,1,""],b:[33,2,1,""],bounds:[33,2,1,""],c:[33,2,1,""],connecting_span:[33,2,1,""],consecutive:[33,2,1,""],context_length:[33,2,1,""],cyclic:[33,2,1,""],differences:[33,2,1,""],flip:[33,2,1,""],force_context:[33,2,1,""],get_slice:[33,2,1,""],get_slice_iter:[33,2,1,""],i:[33,2,1,""],index:[33,2,1,""],intersection:[33,2,1,""],invert:[33,2,1,""],overlaps_with:[33,2,1,""],ranges:[33,2,1,""],reindex:[33,2,1,""],same_context:[33,2,1,""],slices:[33,2,1,""],sub:[33,2,1,""],t:[33,2,1,""]},"dasi.utils.region.Span":{"new":[33,2,1,""],a:[33,2,1,""],b:[33,2,1,""],bounds:[33,2,1,""],c:[33,2,1,""],connecting_span:[33,2,1,""],consecutive:[33,2,1,""],context_length:[33,2,1,""],cyclic:[33,2,1,""],differences:[33,2,1,""],force_context:[33,2,1,""],get_slice:[33,2,1,""],get_slice_iter:[33,2,1,""],i:[33,2,1,""],index:[33,2,1,""],intersection:[33,2,1,""],invert:[33,2,1,""],overlaps_with:[33,2,1,""],ranges:[33,2,1,""],reindex:[33,2,1,""],same_context:[33,2,1,""],slices:[33,2,1,""],sub:[33,2,1,""],t:[33,2,1,""]},"dasi.utils.region.SpanError":{with_traceback:[33,2,1,""]},"dasi.utils.sequence_design":{design_primers:[34,4,1,""],get_primer_extensions:[34,4,1,""]},dasi:{command_line:[2,0,0,"-"],constants:[3,0,0,"-"],cost:[4,0,0,"-"],design:[5,0,0,"-"],exceptions:[8,0,0,"-"],models:[37,0,0,"-"],utils:[40,0,0,"-"]}},objnames:{"0":["py","module","Python module"],"1":["py","class","Python class"],"2":["py","method","Python method"],"3":["py","attribute","Python attribute"],"4":["py","function","Python function"],"5":["py","exception","Python exception"]},objtypes:{"0":"py:module","1":"py:class","2":"py:method","3":"py:attribute","4":"py:function","5":"py:exception"},terms:{"01t15":7,"13t10":7,"abstract":4,"byte":[2,32],"class":[0,2,3,4,5,10,12,14,15,16,17,18,19,20,21,22,23,24,25,26,27,32,33],"default":[2,4,32,33],"float":[4,12],"function":[0,1,4,5,10,11,13,27,28,30,31,32,33,34,40],"import":[36,38],"int":[0,1,2,4,5,12,15,26,27,32,33,34,40],"new":[2,15,16,17,18,20,22,23,24,26,27,32,33],"null":32,"return":[0,1,2,4,5,8,10,11,12,13,15,16,17,18,20,23,24,26,27,29,30,31,32,33,34,40],"static":[0,5,15,27],"true":[4,10,15,19,27,32,33,34,38],"while":33,And:[24,26],For:[10,32,33],Not:[3,32],The:[4,26,32,33,36,38],Then:33,There:33,These:[0,5],Use:[32,33],Used:12,__init__:[15,16,17,18,19,20,21,22,23,24,25],__miss:3,__suffix:32,__traceback__:[8,29,32,33],_alignment_hash:[15,27],_asdict:20,_assemble_graphs_without_thread:[0,5],_batch_add_edge_cost:12,_check_path:13,_create_pcr_product_align:[15,27],_fill_valu:32,_get_iter_non_repeat:[0,5],_get_span_cost:2,_head:19,_make:20,_multinode_to_shortest_path:13,_nodes_to_fullpath:13,_optimize_without_thread:[0,5],_primer_min_span:38,_replac:20,_result:[0,5],_seqdb:[0,5],_share_query_blast:[0,5],_synthesis_left_span_rang:38,_synthesis_step_s:38,a_axi:10,abc:[0,4,5,15,19,26,27,32,33],abs:33,abs_wrap:33,absolut:33,access:33,accord:[15,27,31,32],accumul:30,accur:[15,19,20,21,25],across:[0,4,5,32],activ:25,actual:3,add:[0,5,12,32,33],add_assembl:[0,5],add_edg:12,add_frag:[0,5],add_nod:12,add_prim:[0,5],add_queri:[0,5],add_templ:[0,5],added:32,addit:[12,32],additionalitem:38,additionalproperti:38,adjust:33,aggrag:32,aggreg:32,algorithm:[1,40],alia:20,align:[0,3,5,8,12,15,16,17,18,21,23,24,27],alignment_contain:12,alignment_group:[15,21,27],alignment_hash:[15,27],alignment_typ:[15,27],alignmentcontain:[0,5,11,16,27],alignmentcontainerexcept:8,alignmentcontainerfactori:27,alignmentexcept:8,alignmentgroup:[15,26,27],alignmentgroupbas:[17,23,24,26],all:[0,5,10,12,15,26,27,32,33],allow:[15,27,32],allow_invalid:[0,5],almost:3,along:10,alreadi:3,also:[32,33],alwai:[26,33],amount:38,ani:[1,4,15,27,32,33,34,40],anneal:[4,38],anoth:[32,33],append:32,appli:32,applic:32,applui:32,apply_to_col_nam:32,appropri:34,arang:32,arbitrari:30,arg:[32,33],argument:[10,32],around:[33,34],arr:[13,31,32,33],arrai:[1,10,31,32,38,40],as_typ:33,assembl:[0,4,5,8,11,12,34,37,38],assemble_graph:11,assemblygraph:12,assemblygraphbuild:12,assemblynod:[0,5,12,34],assert:[10,33],asset:33,assign:21,associ:[17,18,23,24,26],assum:[1,40],astyp:32,attribut:[0,3,4,5,15,16,17,18,19,20,22,23,24,25,26,27,32,33],atyp:[12,15,16,17,18,23,24,26,27],automat:33,avail:33,axes:10,axi:[10,32],b_axi:10,back:33,base:[0,2,3,4,5,8,12,15,16,17,18,19,20,21,22,23,24,25,26,27,29,32,33,38],basic:33,batch:12,been:[4,15,27],befor:[4,32],being:27,best:33,between:[12,26,32,33],beyond:4,bind:3,bisect:[1,40],bisect_between:[1,40],bisect_slice_between:[1,40],blast:[1,16,27,40],blast_to_region:27,block:32,bool:[4,33,34],both:[1,26,40],bound:33,broadcast:[10,32],build:[4,11,12],builder:4,built:33,calcul:[4,12,33],call:2,callabl:[1,32,33,40],can:[3,4,10,15,16,24,26,27,32,33],certain:33,chang:[32,33],check:[1,13,40],choos:31,circular:33,classmethod:[0,4,5,15,20,27,32],code:[32,36],col:[10,32],collect:[0,5,15,19,26,27,32,33],column:[10,32,38],combin:[15,27],come:[35,39],command:36,command_lin:36,compil:[0,5],compile_librari:[0,5],complex:32,comput:[4,30],conatin:36,concat:32,concaten:32,condit:[12,31],connect:33,connecting_span:33,consecut:33,consid:[4,38],consider:33,constant:36,construct:[16,27,33],consum:[1,40],contain:[0,5,8,11,12,15,16,19,27,33,36],container_factori:11,container_list:[0,5],content:32,context:33,context_length:33,contribut:36,convers:32,convert:[4,27,33],copi:[10,26,32],cost:[0,2,5,8,12,22,31,36],cost_model:[4,38],costbuild:4,costmodelbas:4,count:20,creat:[15,16,27,32,33],current:32,custom:[4,33,38],cutoff:30,cycl:[1,13,30,40],cycle_endpoint:13,cyclic:[11,13,30,31,33,34],cyclic_sort_kei:30,dai:[4,38],dasi:[36,38],dasicli:2,dasicostparametervalidationerror:8,dasidesignexcept:8,dasiexcept:8,dasiinvalidmolecularassembl:8,dasinoprimerpairsexcept:8,dasisequencedesignexcept:8,dasiwarn:8,data:[1,10,12,16,24,26,27,32,38,40],databas:[16,27],datafram:[4,10,32],decod:4,def:32,defin:10,del:32,delet:32,design:[2,4,8,22,34,36],design_prim:34,designresult:[0,5],dest:12,determin:[1,24,26,40],df1:32,df2:32,df3:32,df_to_np_rang:10,dfs:32,dict:[0,5,15,16,17,18,24,26,27,32,34],dictionari:[15,16,27,34],diff:33,differ:[32,33],digraph:[0,5,11,13,30,34],dimension:32,direct:33,directli:33,directori:2,disallow:[15,27,33],distanc:[12,30],div:32,divid:[31,32],dna:[0,2,3,5,26],do_rais:19,docsrc:2,document:2,doe:[1,40],dtype:[10,32],dump:32,duplic:[10,13],dure:33,each:[0,5,32],edf:4,edg:[12,34],eff_df:4,effect:38,effici:[4,12,22,31,38],either:[4,31],element:[31,32],els:[10,32],emptyspan:33,encod:4,end:[1,15,17,24,26,27,33,34,40],endpoint:[13,33],entir:[1,40],equival:[32,33],error:[29,33],estim:13,evalu:4,event:38,everyth:33,exact:[24,26],exampl:[4,10,32,33,38],except:[32,33,36],exclus:[1,27,33,34,40],exist:[3,4,15,27],expand:[10,15,20,27,32],expand_overlap:[15,27],expand_prim:[15,27],expand_primer_dim:[15,27],expand_primer_pair:[15,27],expans:10,expect:32,explain:34,explicitli:31,ext:4,extend:[3,4],extens:[4,34],extra:[24,26],fals:[0,4,5,13,15,22,27,30,32,33,38],fasta:2,faster:31,field:20,file:8,fill:32,fill_valu:32,filter_linear_record:[0,5],filter_perfect_subject:[0,5],find:[0,1,5,10,30,31,33,40],find_all_min_path:31,find_groups_by_po:[15,27],first:[19,24,26,33],flank:34,flatten:10,flatten_axi:10,flexibl:4,flip:33,float64:10,follow:[4,15,27,32,33],forc:32,force_context:33,form:[15,27],format:[16,27],forward:[23,24,26],found:[4,33],fragment:[0,2,3,5,15,27,34],fragment_shared_with_other_queri:3,frame:[10,32,38],freez:[15,27],from:[0,3,4,5,12,13,15,16,20,24,26,27,31,32,33,36,38],from_json:4,frozendict:[16,27],full:[13,33],full_assembly_graph:19,func:32,further:[15,27],fwd:[3,15,23,24,26,27],gap:3,gener:[3,8,29,32,33],get:[19,32,34],get_groups_by_typ:[15,27],get_primer_extens:34,get_slic:33,get_slice_it:33,gibson:38,github:2,given:[32,34],global:38,global_material_modifi:38,global_time_cost:38,goal:[0,2,5,20],graph:[0,5,11,12,13,19,28,30,34],greater:33,group:[3,15,16,17,18,23,24,26,27,32],group_appli:32,group_typ:[15,17,18,23,24,26,27],groups_by_typ:[15,27],guidelin:36,had:4,halfwai:32,handl:34,has:[10,32,33],hashabl:[15,27],have:[4,15,24,26,27,33],help:[15,19,20,21,25],high:[1,36,40],homolog:38,horizont:32,how:[4,36],howev:32,hstack:[10,32],ident:26,identifi:26,idtprim:38,idtultram:38,ignor:33,ignore_invalid:[0,5],ignore_wrap:33,immedi:33,implement:31,incid:13,includ:10,inclus:[1,26,27,33,40],index:[20,27,30,32,33],index_slic:13,indexerror:33,indic:[1,4,13,27,33,40],indici:[4,17,18,23,24,26,33],individu:32,inersect:33,inf:[3,4,31,38],infer_typ:33,infin:3,inform:[19,36],init:30,initi:[2,4,7,15,16,19,20,21,22,23,24,25,26,27,32,33,38],inplac:32,input:[4,8,25,36],input_data:4,insid:[24,26],instanc:[15,21,27],integ:4,inter:33,interfac:2,intern:34,internal_or_extern:12,intersect:[24,26,33],invalid:[4,8],invers:33,invert:33,item:[13,32,38],iter:[0,1,5,19,20,32,33,40],itertool:33,its:34,job:[0,2,5],json:[4,8,16,27,38],justinvrana:2,kei:[0,1,5,15,16,17,18,23,24,26,27,30,31,32,40],keyword:32,kwarg:[12,32],kwd:20,label:[15,27],last:[24,26,33],left:[3,4,34],left_ext:34,left_overhang:34,left_span_rang:4,len:33,length:[26,33],level:36,librari:[0,5],librarydesign:[0,5],like:[32,33],lim:31,lim_siz:[15,27],line:36,linear:[0,5,33],list:[0,1,4,5,11,13,15,16,17,18,23,24,26,27,30,32,33,40],list_of_arr:32,load:[4,32],load_blast_json:[16,27],locat:[20,32,34],loggabl:[15,27],logger:[15,27],look:[4,34],low:[1,40],lowest:[30,33],lprimer_left_ext:34,lseq:34,made:3,mai:26,maintain:[0,5,16,27,31,32],make:[20,26],mani:32,manipul:34,map:[20,32,33],mata:31,matb:31,materi:[0,4,5,12,38],material_mod:4,material_modifi:4,mathemat:32,matrici:31,max:[10,38],max_col:10,max_homolog:3,max_siz:22,maximum:[3,38],maxitem:38,maxproperti:38,mean:32,memori:32,merg:32,meta:[17,18,24,26],metadata:21,metatyp:22,method:[0,2,4,5,12,15,16,17,18,19,20,21,22,23,24,25,26,27,32,33,34],min:[10,38],min_ann:4,min_col:10,min_overlap:3,min_siz:22,min_span:4,minimum:[3,31,38],minitem:38,minproperti:38,miss:[3,32],mixin:26,model:[36,38],modif:[15,27],modifi:38,modul:[0,2,3,4,5,8,28],molecul:[3,22,25,37],molecular:[4,8,26,37],molecule_typ:21,msgpack:[4,32],multi:32,multipcrproductalignmentgroup:26,multipl:[32,33],multipleof:38,multipli:[4,32,38],multipoint_shortest_path:30,multiprocess:[0,5,11],multiprocessing_assemble_graph:11,multiprocessing_optimize_graph:11,must:[3,26,32],my_paramet:4,n_job:[0,2,5,11],n_path:[0,5,11,13],name:[12,17,18,20,22,24,25,26,32,33,38],nan:[31,32],nativ:32,ndarrai:[4,10,31,32],networkxutilsexcept:29,node:[0,5,12,13,19,30,34],non:[24,26,33],none:[0,4,5,10,12,13,15,16,17,18,21,22,24,25,26,27,30,31,32,33,34],now:33,number:[0,2,4,5,15,20,27,33,38],numpi:[4,10,32],numpydatafram:[4,32],numpydataframeexcept:32,numpydataframeindex:32,obj:4,object:[0,2,3,5,12,16,20,21,22,25,26,27,32,38],occurr:20,onc:12,one:[24,26,33],onli:[0,5,32],onto:33,open:[4,38],oper:[32,33],optim:[0,5,11],optimize_librari:[0,5],option:[4,12,15,17,18,24,26,27,33,34],order:31,ordereddict:20,origin:34,other:[3,15,25,27,32,33,34],ought:[24,26],output:25,outsid:33,overal:12,overhang:[20,24,26,34],overlap:[3,15,27,33],overlaps_with:33,packag:2,pair:[15,26,27,34],panda:[4,10,32],param:[4,32],paramet:[0,1,2,4,5,8,12,13,15,16,17,18,23,24,26,27,30,32,33,34,36,40],parrallel:2,path:[0,1,5,12,13,19,30,31,32,40],pcr:[3,15,23,24,26,27],pcr_product:[3,15,27],pcr_product_with_left_prim:3,pcr_product_with_prim:3,pcr_product_with_right_prim:3,pcrproductalignmentgroup:[15,26,27],pdf:4,per:[4,38],perfect_subject:[1,40],perform:[10,32],plot:4,point:[15,24,26,27,33,34],pos:33,posit:[15,17,26,27,33],possibl:[15,27,33],potenti:3,pre:[1,3,40],predecessor:34,prefix:32,prefix_:32,preprocess:32,primer:[0,2,3,4,5,15,23,24,26,27,34,38],primer_cost:[4,38],primer_cost_model:4,primer_df:4,primer_effici:38,primer_extension_product:3,primer_extension_product_with_left_prim:3,primer_extension_product_with_right_prim:3,primer_min_ann:38,primer_min_bind:3,primer_param:4,primercostmodel:4,primerparam:4,print:[2,32,33],procedur:33,prod:31,produc:[0,3,4,5,15,17,18,19,23,24,25,26,27],product:[3,15,23,24,26,27],project:7,properli:32,properti:[0,4,5,15,16,17,18,20,23,24,26,27,32,33,38],provid:[0,1,2,3,4,5,8,10,24,26,28,32,33,34,37,38,40],qend:[17,18,23,24,26],qstart:[17,18,23,24,26],queri:[0,3,5,11,15,16,17,18,19,23,24,26,27],query_kei:[0,5,16,17,18,19,23,24,26,27],query_length:11,query_or_subject:27,query_region:[17,18,21,23,24,26],rais:[32,33],rang:[4,33],reaction:[24,26],record:[0,5,16,26,27],recov:13,reduc:33,reduce_op:33,redund:[23,26],redundent_alignment_group:[15,27],reflect:33,region:[0,3,5,15,17,18,23,24,26,27,34],region_id:33,registri:[0,5],reindex:33,rel:26,relev:19,reload:2,remain:33,remaining_col:10,remap:33,repeat:[0,5],replac:[20,31],repr:32,repres:[3,4,15,17,18,20,24,26,27,32,33],represent:[15,27,37],representsmolecul:[18,26],reshap:32,restor:34,result:[0,1,5,16,27,32,33,34,40],retriev:[16,27],return_length:13,reus:3,rev:[3,15,23,24,26,27],revers:[23,24,26],right:[3,4,34],right_ext:34,right_overhang:34,row:32,rprimer_right_ext:34,rseq:34,rtype:32,run:[2,32],same:[15,17,26,27,31,32,33],same_context:33,save:2,schema:4,sdf:4,search:[3,15,27],see:[15,19,20,21,25],select:[31,32],select_from_arr:31,self:[8,15,19,20,21,25,29,30,32,33],seqdb:[0,5,15,16,19,27],seqrecord:[0,5,15,16,27],sequenc:[0,5,8,16,17,20,21,26,27,34],set:[0,4,5,8,15,27,29,32,33],sever:[24,25,26,32,33],shallow:26,shape:[10,31,32],share:[0,3,5,15,16,17,26,27],shared_frag:3,shortest:[13,30],should:19,shown:33,side:34,signatur:[15,19,20,21,25],similarli:33,simpli:33,singl:[15,27,32],situat:[15,24,26,27],size:[4,15,26,27,33],slice:[1,32,33,40],somehow:31,soon:[35,39],sort:[1,30,31,40],sort_cycl:31,sort_with_kei:[1,40],sourc:[0,1,2,3,4,5,8,10,11,12,13,15,16,17,18,19,21,22,23,24,25,26,27,29,30,31,32,33,34,40],span:[0,2,4,5,10,12,33,34],span_cost:[0,5,11,12],spancost:[0,4,5,38],spanerror:33,special:[1,40],specif:10,specifi:[1,10,11,13,17,18,20,23,24,26,32,34,40],square_broadcast:10,src:12,stack:32,start:[1,15,17,26,27,33,40],startpoint:33,step:4,step_siz:4,still:33,str:[0,5,12,15,16,17,18,23,24,26,27,30,32,33,34],strategi:32,strict:33,string:[32,34,38],sub:[17,18,23,24,26,33],sub_region:[17,18,23,24,26],subject:[0,1,5,15,17,18,23,24,26,27,40],subject_kei:[17,18,23,24,26],subject_region:[17,18,23,24,26],subregion:26,successor:34,suffix:32,sum:31,summari:33,summary_func:33,support:30,sympy_dijkstra:30,syn_cost:4,synthes:[3,22],synthesi:4,synthesis_cost:38,synthesis_df:4,synthesiscostmodel:4,synthesisparam:4,t_co:[1,40],take:[4,19,25,33],taken:26,target:30,templat:[0,2,3,5,23,24,26,34],template_group:[15,27],terrarium:29,than:[31,33],thei:32,them:32,themselv:[0,5],thi:[0,1,2,3,4,5,8,15,17,18,23,24,26,27,28,31,32,33,34,37,38,40],those:[15,27],three:[15,27],through:[30,33],throw_error:33,thrown:33,time:[4,12,33,38],time_cost:4,to_df:[4,32],togeth:[15,27,32],top:[0,5],topolog:33,touch:13,translat:33,treat:33,tupl:[0,1,4,5,11,13,15,20,27,32,33,34,40],two:[3,12,26,27,31,32,33],type:[0,1,4,5,11,12,13,15,16,17,18,19,20,21,23,24,25,26,27,32,33,34,38,40],typecast:33,underli:[32,33],unfreez:[15,27],unintuit:33,union:[4,15,24,26,27,33,34],uniqueitem:38,unlik:32,unmap:33,updat:32,usag:[32,36],use:[0,1,4,5,30,40],use_direct:22,used:[10,32,33],user:2,uses:[3,33],using:[1,11,15,27,30,32,33,38,40],util:36,valid:[4,13,15,26,27,32,33,38],valu:[20,32],variou:[1,28,40],verbos:2,version:2,via:[16,27],view:32,vstack:32,wai:33,wait:38,warn:8,weight:[30,31],weight_kei:30,well:33,whatev:31,when:[2,33],whenev:33,whether:[1,4,30,33,34,40],which:[4,10,20,21,24,26],whose:[0,5],wise:32,with_traceback:[8,29,32,33],within:[10,15,24,26,27,33],would:[4,32,33],wrap:33,you:32,zip:13},titles:["dasi.design","dasi.utils","Command Line (<code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">dasi.command_line</span></code>)","Constants (<code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">dasi.constants</span></code>)","Cost Model (<code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">dasi.cost</span></code>)","Design (<code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">dasi.design</span></code>)","&lt;no title&gt;","dasi change log","Exceptions (<code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">dasi.exceptions</span></code>)","dasi.cost.params","dasi.cost.utils","dasi.design.design_algorithms","dasi.design.graph_builder","dasi.design.optimize","dasi.design.plotter","dasi.models.AlignmentContainer","dasi.models.AlignmentContainerFactory","dasi.models.AlignmentGroup","dasi.models.AlignmentGroupBase","dasi.models.Assembly","dasi.models.AssemblyNode","dasi.models.Molecule","dasi.models.MoleculeType","dasi.models.MultiPCRProductAlignmentGroup","dasi.models.PCRProductAlignmentGroup","dasi.models.Reaction","dasi.models.alignment","dasi.models.alignment_container","dasi.utils.networkx","dasi.utils.networkx.exceptions","dasi.utils.networkx.shortest_path","dasi.utils.networkx.utils","dasi.utils.npdf","dasi.utils.region","dasi.utils.sequence_design","Code Guidelines","API Reference","Models (<code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">dasi.models</span></code>)","User Parameter Inputs","Usage","Utilities (<code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">dasi.utils</span></code>)"],titleterms:{"default":38,align:26,alignment_contain:27,alignmentcontain:15,alignmentcontainerfactori:16,alignmentgroup:17,alignmentgroupbas:18,api:36,assembl:19,assemblynod:20,chang:7,code:35,command:2,command_lin:2,constant:3,cost:[4,9,10,38],dasi:[0,1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,37,40],design:[0,5,11,12,13,14],design_algorithm:11,develop:36,document:36,except:[8,29],graph_build:12,guidelin:35,input:38,line:2,log:7,model:[4,15,16,17,18,19,20,21,22,23,24,25,26,27,37],modul:[1,37,40],molecul:21,moleculetyp:22,multipcrproductalignmentgroup:23,networkx:[1,28,29,30,31,40],npdf:32,optim:13,param:9,paramet:38,pcrproductalignmentgroup:24,plotter:14,reaction:25,refer:36,region:33,schema:38,sequence_design:34,shortest_path:30,usag:39,user:[36,38],util:[1,4,10,28,29,30,31,32,33,34,40]}})