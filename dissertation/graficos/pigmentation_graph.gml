graph [
  directed 1
  multigraph 1
  typedefs [
    id "ends_during"
    name "ends during"
    namespace "external"
    xref "_networkx_list_start"
    xref "RO:0002093"
  ]
  typedefs [
    id "happens_during"
    name "happens during"
    namespace "external"
    xref "_networkx_list_start"
    xref "RO:0002092"
    is_transitive "true"
    is_a "_networkx_list_start"
    is_a "ends_during"
  ]
  typedefs [
    id "has_ontology_root_term"
    name "has ontology root term"
    namespace "external"
    xref "_networkx_list_start"
    xref "IAO:0000700"
    is_metadata_tag "true"
    is_class_level "true"
  ]
  typedefs [
    id "has_part"
    name "has part"
    namespace "external"
    xref "_networkx_list_start"
    xref "BFO:0000051"
    is_transitive "true"
  ]
  typedefs [
    id "negatively_regulates"
    name "negatively regulates"
    namespace "external"
    xref "_networkx_list_start"
    xref "RO:0002212"
    is_a "_networkx_list_start"
    is_a "regulates"
  ]
  typedefs [
    id "occurs_in"
    name "occurs in"
    namespace "external"
    xref "_networkx_list_start"
    xref "BFO:0000066"
    transitive_over "_networkx_list_start"
    transitive_over "part_of"
  ]
  typedefs [
    id "part_of"
    name "part of"
    namespace "external"
    xref "_networkx_list_start"
    xref "BFO:0000050"
    is_transitive "true"
    inverse_of "_networkx_list_start"
    inverse_of "has_part"
  ]
  typedefs [
    id "positively_regulates"
    name "positively regulates"
    namespace "external"
    xref "_networkx_list_start"
    xref "RO:0002213"
    holds_over_chain "_networkx_list_start"
    holds_over_chain "negatively_regulates negatively_regulates"
    is_a "_networkx_list_start"
    is_a "regulates"
  ]
  typedefs [
    id "regulates"
    name "regulates"
    namespace "external"
    xref "_networkx_list_start"
    xref "RO:0002211"
    is_transitive "true"
  ]
  typedefs [
    id "starts_during"
    name "starts during"
    namespace "external"
    xref "_networkx_list_start"
    xref "RO:0002091"
  ]
  typedefs [
    id "term_tracker_item"
    name "term tracker item"
    namespace "external"
    xref "_networkx_list_start"
    xref "IAO:0000233"
    is_metadata_tag "true"
    is_class_level "true"
  ]
  instances "[]"
  node [
    id 0
    label "GO:0033059"
    name "cellular pigmentation"
    namespace "biological_process"
    def "&#34;The deposition or aggregation of coloring matter in a cell.&#34; [GOC:mtg_MIT_16mar07]"
    is_a "GO:0009987"
    is_a "GO:0043473"
    enriched 1
    fdr 0.001106195948408
  ]
  node [
    id 1
    label "GO:0043473"
    name "pigmentation"
    namespace "biological_process"
    def "&#34;The accumulation of pigment in an organism, tissue or cell, either by increased deposition or by increased number of cells.&#34; [GOC:jl]"
    subset "goslim_chembl"
    subset "goslim_pir"
    is_a "_networkx_list_start"
    is_a "GO:0008150"
  ]
  node [
    id 2
    label "GO:0005737"
    name "cytoplasm"
    namespace "cellular_component"
    def "&#34;The contents of a cell excluding the plasma membrane and nucleus, but including other subcellular structures.&#34; [ISBN:0198547684]"
    subset "goslim_candida"
    subset "goslim_chembl"
    subset "goslim_metagenomics"
    subset "goslim_pir"
    subset "goslim_plant"
    subset "goslim_yeast"
    subset "prokaryote_subset"
    xref "_networkx_list_start"
    xref "Wikipedia:Cytoplasm"
    is_a "_networkx_list_start"
    is_a "GO:0110165"
    relationship "_networkx_list_start"
    relationship "part_of GO:0005622"
    property_value "_networkx_list_start"
    property_value "term_tracker_item https://github.com/geneontology/go-ontology/issues/23023 xsd:anyURI"
  ]
  node [
    id 3
    label "GO:0005575"
    name "cellular_component"
    namespace "cellular_component"
    alt_id "_networkx_list_start"
    alt_id "GO:0008372"
    def "&#34;A location, relative to cellular compartments and structures, occupied by a macromolecular machine when it carries out a molecular function. There are two ways in which the gene ontology describes locations of gene products: (1) relative to cellular structures (e.g., cytoplasmic side of plasma membrane) or compartments (e.g., mitochondrion), and (2) the stable macromolecular complexes of which they are parts (e.g., the ribosome).&#34; [GOC:pdt, NIF_Subcellular:sao1337158144]"
    comment "Note that, in addition to forming the root of the cellular component ontology, this term is recommended for use for the annotation of gene products whose cellular component is unknown. When this term is used for annotation, it indicates that no information was available about the cellular component of the gene product annotated as of the date the annotation was made; the evidence code 'no data' (ND), is used to indicate this."
    subset "goslim_candida"
    subset "goslim_chembl"
    subset "goslim_metagenomics"
    subset "goslim_pir"
    subset "goslim_plant"
    subset "goslim_yeast"
    synonym "&#34;cell or subcellular entity&#34; EXACT []"
    synonym "&#34;cellular component&#34; EXACT []"
    synonym "&#34;subcellular entity&#34; RELATED [NIF_Subcellular:nlx_subcell_100315]"
    xref "_networkx_list_start"
    xref "NIF_Subcellular:sao1337158144"
  ]
  node [
    id 4
    label "GO:0008150"
    name "biological_process"
    namespace "biological_process"
    alt_id "GO:0000004"
    alt_id "GO:0007582"
    alt_id "GO:0044699"
    def "&#34;A biological process is the execution of a genetically-encoded biological module or program. It consists of all the steps required to achieve the specific biological objective of the module. A biological process is accomplished by a particular set of molecular functions carried out by specific gene products (or macromolecular complexes), often in a highly regulated manner and in a particular temporal sequence.&#34; [GOC:pdt]"
    comment "Note that, in addition to forming the root of the biological process ontology, this term is recommended for use for the annotation of gene products whose biological process is unknown. When this term is used for annotation, it indicates that no information was available about the biological process of the gene product annotated as of the date the annotation was made; the evidence code 'no data' (ND), is used to indicate this."
    subset "goslim_candida"
    subset "goslim_chembl"
    subset "goslim_metagenomics"
    subset "goslim_pir"
    subset "goslim_plant"
    subset "goslim_pombe"
    subset "goslim_yeast"
    synonym "&#34;biological process&#34; EXACT []"
    synonym "&#34;physiological process&#34; EXACT []"
    synonym "&#34;single organism process&#34; RELATED []"
    synonym "&#34;single-organism process&#34; RELATED []"
    xref "_networkx_list_start"
    xref "Wikipedia:Biological_process"
    property_value "_networkx_list_start"
    property_value "term_tracker_item https://github.com/geneontology/go-ontology/issues/24968 xsd:anyURI"
    created_by "jl"
    creation_date "2012-09-19T15:05:24Z"
  ]
  node [
    id 5
    label "GO:0043226"
    name "organelle"
    namespace "cellular_component"
    def "&#34;Organized structure of distinctive morphology and function. Includes the nucleus, mitochondria, plastids, vacuoles, vesicles, ribosomes and the cytoskeleton, and prokaryotic structures such as anammoxosomes and pirellulosomes. Excludes the plasma membrane.&#34; [GOC:go_curators]"
    subset "goslim_chembl"
    subset "goslim_generic"
    subset "goslim_pir"
    subset "prokaryote_subset"
    xref "NIF_Subcellular:sao1539965131"
    xref "Wikipedia:Organelle"
    is_a "_networkx_list_start"
    is_a "GO:0110165"
  ]
  node [
    id 6
    label "GO:0031410"
    name "cytoplasmic vesicle"
    namespace "cellular_component"
    alt_id "_networkx_list_start"
    alt_id "GO:0016023"
    def "&#34;A vesicle found in the cytoplasm of a cell.&#34; [GOC:ai, GOC:mah, GOC:vesicles]"
    subset "goslim_agr"
    subset "goslim_candida"
    subset "goslim_chembl"
    subset "goslim_generic"
    subset "goslim_mouse"
    subset "goslim_yeast"
    synonym "&#34;cytoplasmic membrane bounded vesicle&#34; RELATED []"
    synonym "&#34;cytoplasmic membrane-enclosed vesicle&#34; RELATED []"
    synonym "&#34;cytoplasmic, membrane-bounded vesicle&#34; RELATED []"
    xref "_networkx_list_start"
    xref "NIF_Subcellular:sao180601769"
    is_a "_networkx_list_start"
    is_a "GO:0097708"
    intersection_of "GO:0031982"
    intersection_of "part_of GO:0005737"
    relationship "_networkx_list_start"
    relationship "part_of GO:0005737"
  ]
  node [
    id 7
    label "GO:0048770"
    name "pigment granule"
    namespace "cellular_component"
    def "&#34;A small, subcellular membrane-bounded vesicle containing pigment and/or pigment precursor molecules. Pigment granule biogenesis is poorly understood, as pigment granules are derived from multiple sources including the endoplasmic reticulum, coated vesicles, lysosomes, and endosomes.&#34; [GOC:jid, GOC:mh]"
    is_a "_networkx_list_start"
    is_a "GO:0031410"
    enriched 1
    fdr 0.0002552980472753
  ]
  node [
    id 8
    label "GO:0031982"
    name "vesicle"
    namespace "cellular_component"
    alt_id "_networkx_list_start"
    alt_id "GO:0031988"
    def "&#34;Any small, fluid-filled, spherical organelle enclosed by membrane.&#34; [GOC:mah, GOC:pz, GOC:vesicles]"
    subset "_networkx_list_start"
    subset "goslim_pir"
    synonym "&#34;membrane-bounded vesicle&#34; RELATED []"
    synonym "&#34;membrane-enclosed vesicle&#34; RELATED []"
    xref "NIF_Subcellular:sao221389602"
    xref "Wikipedia:Vesicle_(biology)"
    is_a "_networkx_list_start"
    is_a "GO:0043227"
  ]
  node [
    id 9
    label "GO:0005622"
    name "intracellular anatomical structure"
    namespace "cellular_component"
    def "&#34;A component of a cell contained within (but not including) the plasma membrane. In eukaryotes it includes the nucleus and cytoplasm.&#34; [ISBN:0198506732]"
    subset "gocheck_do_not_annotate"
    subset "goslim_chembl"
    subset "goslim_metagenomics"
    subset "goslim_plant"
    synonym "&#34;internal to cell&#34; EXACT []"
    synonym "&#34;intracellular&#34; EXACT []"
    synonym "&#34;nucleocytoplasm&#34; RELATED [GOC:mah]"
    synonym "&#34;protoplasm&#34; EXACT []"
    synonym "&#34;protoplast&#34; RELATED [GOC:mah]"
    xref "_networkx_list_start"
    xref "Wikipedia:Intracellular"
    is_a "_networkx_list_start"
    is_a "GO:0110165"
    property_value "_networkx_list_start"
    property_value "term_tracker_item https://github.com/geneontology/go-ontology/issues/17776 xsd:anyURI"
  ]
  node [
    id 10
    label "GO:0043229"
    name "intracellular organelle"
    namespace "cellular_component"
    def "&#34;Organized structure of distinctive morphology and function, occurring within the cell. Includes the nucleus, mitochondria, plastids, vacuoles, vesicles, ribosomes and the cytoskeleton. Excludes the plasma membrane.&#34; [GOC:go_curators]"
    subset "_networkx_list_start"
    subset "goslim_pir"
    is_a "_networkx_list_start"
    is_a "GO:0043226"
    intersection_of "GO:0043226"
    intersection_of "part_of GO:0005622"
    relationship "_networkx_list_start"
    relationship "part_of GO:0005622"
  ]
  node [
    id 11
    label "GO:0110165"
    name "cellular anatomical entity"
    namespace "cellular_component"
    def "&#34;A part of a cellular organism that is either an immaterial entity or a material entity with granularity above the level of a protein complex but below that of an anatomical system. Or, a substance produced by a cellular organism with granularity above the level of a protein complex.&#34; [GOC:kmv]"
    subset "gocheck_do_not_annotate"
    subset "goslim_pir"
    is_a "_networkx_list_start"
    is_a "GO:0005575"
    property_value "_networkx_list_start"
    property_value "term_tracker_item https://github.com/geneontology/go-ontology/issues/24200 xsd:anyURI"
    created_by "kmv"
    creation_date "2019-08-12T18:01:37Z"
  ]
  node [
    id 12
    label "GO:0042470"
    name "melanosome"
    namespace "cellular_component"
    def "&#34;A tissue-specific, membrane-bounded cytoplasmic organelle within which melanin pigments are synthesized and stored. Melanosomes are synthesized in melanocyte cells.&#34; [GOC:jl, PMID:11584301]"
    xref "_networkx_list_start"
    xref "Wikipedia:Melanosome"
    is_a "_networkx_list_start"
    is_a "GO:0048770"
    property_value "RO:0002161 NCBITaxon:4895"
    property_value "RO:0002161 NCBITaxon:6237"
    enriched 1
    fdr 0.0002552980472753
  ]
  node [
    id 13
    label "GO:0097708"
    name "intracellular vesicle"
    namespace "cellular_component"
    def "&#34;Any vesicle that is part of the intracellular region.&#34; [GOC:vesicles]"
    is_a "GO:0031982"
    is_a "GO:0043231"
    intersection_of "GO:0031982"
    intersection_of "part_of GO:0005622"
    created_by "pr"
    creation_date "2016-03-29T17:39:45Z"
  ]
  node [
    id 14
    label "GO:0009987"
    name "cellular process"
    namespace "biological_process"
    alt_id "GO:0008151"
    alt_id "GO:0044763"
    alt_id "GO:0050875"
    def "&#34;Any process that is carried out at the cellular level, but not necessarily restricted to a single cell. For example, cell communication occurs among more than one cell, but occurs at the cellular level.&#34; [GOC:go_curators, GOC:isa_complete]"
    comment "This term should not be used for direct annotation. It should be possible to make a more specific annotation to one of the children of this term."
    subset "gocheck_do_not_annotate"
    subset "goslim_plant"
    synonym "&#34;cell growth and/or maintenance&#34; NARROW []"
    synonym "&#34;cell physiology&#34; EXACT []"
    synonym "&#34;cellular physiological process&#34; EXACT []"
    synonym "&#34;single-organism cellular process&#34; RELATED []"
    is_a "_networkx_list_start"
    is_a "GO:0008150"
    disjoint_from "_networkx_list_start"
    disjoint_from "GO:0044848"
    created_by "jl"
    creation_date "2012-12-11T16:56:55Z"
  ]
  node [
    id 15
    label "GO:0043227"
    name "membrane-bounded organelle"
    namespace "cellular_component"
    def "&#34;Organized structure of distinctive morphology and function, bounded by a single or double lipid bilayer membrane. Includes the nucleus, mitochondria, plastids, vacuoles, and vesicles. Excludes the plasma membrane.&#34; [GOC:go_curators]"
    synonym "_networkx_list_start"
    synonym "&#34;membrane-enclosed organelle&#34; EXACT []"
    xref "_networkx_list_start"
    xref "NIF_Subcellular:sao414196390"
    is_a "_networkx_list_start"
    is_a "GO:0043226"
    intersection_of "GO:0043226"
    intersection_of "has_part GO:0016020"
  ]
  node [
    id 16
    label "GO:0043231"
    name "intracellular membrane-bounded organelle"
    namespace "cellular_component"
    def "&#34;Organized structure of distinctive morphology and function, bounded by a single or double lipid bilayer membrane and occurring within the cell. Includes the nucleus, mitochondria, plastids, vacuoles, and vesicles. Excludes the plasma membrane.&#34; [GOC:go_curators]"
    subset "_networkx_list_start"
    subset "goslim_pir"
    synonym "_networkx_list_start"
    synonym "&#34;intracellular membrane-enclosed organelle&#34; EXACT []"
    is_a "GO:0043227"
    is_a "GO:0043229"
    intersection_of "GO:0043227"
    intersection_of "part_of GO:0005622"
  ]
  edge [
    source 0
    target 14
    key "is_a"
  ]
  edge [
    source 0
    target 1
    key "is_a"
  ]
  edge [
    source 1
    target 4
    key "is_a"
  ]
  edge [
    source 2
    target 11
    key "is_a"
  ]
  edge [
    source 2
    target 9
    key "part_of"
  ]
  edge [
    source 5
    target 11
    key "is_a"
  ]
  edge [
    source 6
    target 13
    key "is_a"
  ]
  edge [
    source 6
    target 2
    key "part_of"
  ]
  edge [
    source 7
    target 6
    key "is_a"
  ]
  edge [
    source 8
    target 15
    key "is_a"
  ]
  edge [
    source 9
    target 11
    key "is_a"
  ]
  edge [
    source 10
    target 5
    key "is_a"
  ]
  edge [
    source 10
    target 9
    key "part_of"
  ]
  edge [
    source 11
    target 3
    key "is_a"
  ]
  edge [
    source 12
    target 7
    key "is_a"
  ]
  edge [
    source 13
    target 8
    key "is_a"
  ]
  edge [
    source 13
    target 16
    key "is_a"
  ]
  edge [
    source 14
    target 4
    key "is_a"
  ]
  edge [
    source 15
    target 5
    key "is_a"
  ]
  edge [
    source 16
    target 15
    key "is_a"
  ]
  edge [
    source 16
    target 10
    key "is_a"
  ]
]
