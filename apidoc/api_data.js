define({ "api": [
  {
    "type": "get",
    "url": "/mapping",
    "title": "Get a patients sample ID from their clinical ID",
    "name": "GetMapping",
    "group": "Backend",
    "description": "<p>Get a patients sample ID from their clinical ID</p>",
    "query": [
      {
        "group": "Query",
        "type": "String",
        "optional": false,
        "field": "id",
        "description": "<ul> <li>The patient ID</li> </ul>"
      }
    ],
    "header": {
      "fields": {
        "Header": [
          {
            "group": "Header",
            "type": "String",
            "optional": false,
            "field": "Authorization",
            "description": "<p>JWT &quot;Bearer: <token>&quot;</p>"
          }
        ]
      }
    },
    "success": {
      "fields": {
        "Success 200": [
          {
            "group": "Success 200",
            "type": "Object",
            "optional": false,
            "field": "res.mapping",
            "description": "<ul> <li>Patient mapping</li> </ul>"
          }
        ]
      }
    },
    "version": "0.0.0",
    "filename": "backend/search.py",
    "groupTitle": "Backend"
  },
  {
    "type": "get",
    "url": "/patient",
    "title": "Get specific patient and related ones",
    "name": "GetPatient",
    "group": "Backend",
    "description": "<p>Get a specific patient from the clinical DB by ID as well as similar patients to them</p>",
    "query": [
      {
        "group": "Query",
        "type": "String",
        "optional": false,
        "field": "id",
        "description": "<ul> <li>The patient ID to search by</li> </ul>"
      },
      {
        "group": "Query",
        "type": "String",
        "optional": false,
        "field": "fuzz",
        "description": "<ul> <li>The distance from the patient to find similar patients</li> </ul>"
      }
    ],
    "header": {
      "fields": {
        "Header": [
          {
            "group": "Header",
            "type": "String",
            "optional": false,
            "field": "Authorization",
            "description": "<p>JWT &quot;Bearer: <token>&quot;</p>"
          }
        ]
      }
    },
    "success": {
      "fields": {
        "Success 200": [
          {
            "group": "Success 200",
            "type": "Object",
            "optional": false,
            "field": "res.patient",
            "description": "<ul> <li>Patient</li> </ul>"
          }
        ]
      }
    },
    "version": "0.0.0",
    "filename": "backend/search.py",
    "groupTitle": "Backend"
  },
  {
    "type": "get",
    "url": "/graph/like",
    "title": "Get the graph of patient",
    "name": "GetPatientGraph",
    "group": "Backend",
    "description": "<p>Get the nodes and vertices of /patient</p>",
    "query": [
      {
        "group": "Query",
        "type": "String",
        "optional": false,
        "field": "id",
        "description": "<ul> <li>The patient ID to search by</li> </ul>"
      },
      {
        "group": "Query",
        "type": "String",
        "optional": false,
        "field": "fuzz",
        "description": "<ul> <li>The distance from the patient to find similar patients</li> </ul>"
      }
    ],
    "header": {
      "fields": {
        "Header": [
          {
            "group": "Header",
            "type": "String",
            "optional": false,
            "field": "Authorization",
            "description": "<p>JWT &quot;Bearer: <token>&quot;</p>"
          }
        ]
      }
    },
    "success": {
      "fields": {
        "Success 200": [
          {
            "group": "Success 200",
            "type": "Object",
            "optional": false,
            "field": "graph",
            "description": "<ul> <li>Graph</li> </ul>"
          }
        ]
      }
    },
    "version": "0.0.0",
    "filename": "backend/search.py",
    "groupTitle": "Backend"
  },
  {
    "type": "get",
    "url": "/hpo_search",
    "title": "Search for patients by related HPO term",
    "name": "HPOSearch",
    "group": "Backend",
    "description": "<p>Search for patients by related HPO term</p>",
    "query": [
      {
        "group": "Query",
        "type": "String",
        "optional": false,
        "field": "hpoTerm",
        "description": "<ul> <li>The HPO term to search patients by</li> </ul>"
      },
      {
        "group": "Query",
        "type": "String",
        "optional": false,
        "field": "fuzz",
        "description": "<ul> <li>The distance from hpoTerm to search</li> </ul>"
      },
      {
        "group": "Query",
        "type": "String",
        "optional": false,
        "field": "cutoff",
        "description": "<ul> <li>How many patients need to have the phenotype to be included</li> </ul>"
      },
      {
        "group": "Query",
        "type": "String",
        "optional": false,
        "field": "limitation",
        "description": "<ul> <li>A HPO term all patients must have</li> </ul>"
      }
    ],
    "header": {
      "fields": {
        "Header": [
          {
            "group": "Header",
            "type": "String",
            "optional": false,
            "field": "Authorization",
            "description": "<p>JWT &quot;Bearer: <token>&quot;</p>"
          }
        ]
      }
    },
    "success": {
      "fields": {
        "Success 200": [
          {
            "group": "Success 200",
            "type": "Object[]",
            "optional": false,
            "field": "res.overall_results",
            "description": "<ul> <li>All HPO terms found and their patients</li> </ul>"
          },
          {
            "group": "Success 200",
            "type": "Object[]",
            "optional": false,
            "field": "res.patients",
            "description": "<ul> <li>All patients found</li> </ul>"
          },
          {
            "group": "Success 200",
            "type": "Object[]",
            "optional": false,
            "field": "res.results_per_dataset",
            "description": "<ul> <li>All HPO terms found and their patients, group by cohort</li> </ul>"
          },
          {
            "group": "Success 200",
            "type": "Object[]",
            "optional": false,
            "field": "res.cohortsDenied",
            "description": "<ul> <li>Cohorts that had patients, but the JWT did not have access for</li> </ul>"
          }
        ]
      }
    },
    "version": "0.0.0",
    "filename": "backend/search.py",
    "groupTitle": "Backend"
  },
  {
    "type": "get",
    "url": "/graph/with",
    "title": "Get the graph of hpo_search",
    "name": "HPOSearchGraph",
    "group": "Backend",
    "description": "<p>Get the nodes and vertices of /hpo_search</p>",
    "query": [
      {
        "group": "Query",
        "type": "String",
        "optional": false,
        "field": "hpoTerm",
        "description": "<ul> <li>The HPO term to search patients by</li> </ul>"
      },
      {
        "group": "Query",
        "type": "String",
        "optional": false,
        "field": "fuzz",
        "description": "<ul> <li>The distance from hpoTerm to search</li> </ul>"
      },
      {
        "group": "Query",
        "type": "String",
        "optional": false,
        "field": "cutoff",
        "description": "<ul> <li>How many patients need to have the phenotype to be included</li> </ul>"
      },
      {
        "group": "Query",
        "type": "String",
        "optional": false,
        "field": "limitation",
        "description": "<ul> <li>A HPO term all patients must have</li> </ul>"
      }
    ],
    "header": {
      "fields": {
        "Header": [
          {
            "group": "Header",
            "type": "String",
            "optional": false,
            "field": "Authorization",
            "description": "<p>JWT &quot;Bearer: <token>&quot;</p>"
          }
        ]
      }
    },
    "success": {
      "fields": {
        "Success 200": [
          {
            "group": "Success 200",
            "type": "Object",
            "optional": false,
            "field": "graph",
            "description": "<ul> <li>Graph</li> </ul>"
          }
        ]
      }
    },
    "version": "0.0.0",
    "filename": "backend/search.py",
    "groupTitle": "Backend"
  },
  {
    "type": "get",
    "url": "/variant_search",
    "title": "Find patients by variants",
    "name": "VariantSearch",
    "group": "Backend",
    "description": "<p>Search for patients by variants</p>",
    "query": [
      {
        "group": "Query",
        "type": "String",
        "optional": false,
        "field": "region",
        "description": "<ul> <li>CSV of regions</li> </ul>"
      },
      {
        "group": "Query",
        "type": "String",
        "optional": true,
        "field": "ref",
        "description": "<ul> <li>The reference allele of the variant</li> </ul>"
      },
      {
        "group": "Query",
        "type": "String",
        "optional": true,
        "field": "alt",
        "description": "<ul> <li>The alternate allel of the variant</li> </ul>"
      }
    ],
    "header": {
      "fields": {
        "Header": [
          {
            "group": "Header",
            "type": "String",
            "optional": false,
            "field": "Authorization",
            "description": "<p>JWT &quot;Bearer: <token>&quot;</p>"
          }
        ]
      }
    },
    "success": {
      "fields": {
        "Success 200": [
          {
            "group": "Success 200",
            "type": "Object[]",
            "optional": false,
            "field": "res.patient_info",
            "description": "<ul> <li>Patients found</li> </ul>"
          },
          {
            "group": "Success 200",
            "type": "Object[]",
            "optional": false,
            "field": "res.phenotype_count",
            "description": "<ul> <li>Number of phenotypes in those patients found</li> </ul>"
          },
          {
            "group": "Success 200",
            "type": "Object[]",
            "optional": false,
            "field": "res.total_patients",
            "description": "<ul> <li>Number of patients found</li> </ul>"
          },
          {
            "group": "Success 200",
            "type": "Object[]",
            "optional": false,
            "field": "res.variants",
            "description": "<ul> <li>List of variants in patients found in regions</li> </ul>"
          },
          {
            "group": "Success 200",
            "type": "Object[]",
            "optional": false,
            "field": "res.per_cohort",
            "description": "<ul> <li>Data organised per cohort</li> </ul>"
          }
        ]
      }
    },
    "version": "0.0.0",
    "filename": "backend/search.py",
    "groupTitle": "Backend"
  }
] });
