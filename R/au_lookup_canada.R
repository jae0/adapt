#' @title au_lookup_canada
#' @description Lookup tables between StatsCan areal units and Covid19's.
#' @param returntable default is FALSE
#' @param statscan vector of statscan health units, forcing return of matching Covid19 areal units.
#' @param covid19 vector of Covid19 areal units, forcing return of matching statscan health units.
#' @return  Lookup tables between StatsCan areal units and Covid19's.
#' @author Jae Choi, \email{choi.jae.seok@gmail.com}
#' @export
au_lookup_canada = function( returntable=FALSE, statscan=NULL, covid19=NULL) {

  au_look =  as.data.frame( matrix( c(
    "1",	"Newfoundland and Labrador __ Eastern Regional Health Authority",	"NL __ Eastern",
    "2",	"Newfoundland and Labrador __ Central Regional Health Authority",	"NL __ Central",
    "3",	"Newfoundland and Labrador __ Western Regional Health Authority",	"NL __ Western",
    "4",	"Newfoundland and Labrador __ Labrador-Grenfell Regional Health Authority",	"NL __ Labrador-Grenfell",
    "5",	"Nova Scotia __ Zone 1 - Western",	"Nova Scotia __ Zone 1 - Western",
    "6",	"Nova Scotia __ Zone 2 - Northern",	"Nova Scotia __ Zone 2 - Northern",
    "7",	"Nova Scotia __ Zone 3 - Eastern",	"Nova Scotia __ Zone 3 - Eastern",
    "8",	"Nova Scotia __ Zone 4 - Central",	"Nova Scotia __ Zone 4 - Central",
    "9",	"New Brunswick __ Zone 1 (Moncton area)"  ,	"New Brunswick __ Zone 1 (Moncton area)" ,
    "10",	"New Brunswick __ Zone 2 (Saint John area)" ,	"New Brunswick __ Zone 2 (Saint John area)",
    "11",	"New Brunswick __ Zone 3 (Fredericton area)",	"New Brunswick __ Zone 3 (Fredericton area)",
    "12",	"New Brunswick __ Zone 4 (Edmundston area)" ,	"New Brunswick __ Zone 4 (Edmundston area)",
    "13",	"New Brunswick __ Zone 5 (Campbellton area)",	"New Brunswick __ Zone 5 (Campbellton area)",
    "14",	"New Brunswick __ Zone 6 (Bathurst area)" ,	"New Brunswick __ Zone 6 (Bathurst area)",
    "15",	"New Brunswick __ Zone 7 (Miramichi area)" ,	"New Brunswick __ Zone 7 (Miramichi area)",
    "16",	"Quebec __ Région du Bas-Saint-Laurent",	"Quebec __ Bas-Saint-Laurent",
    "17",	"Quebec __ Région du Saguenay - Lac-Saint-Jean",	"Quebec __ Saguenay",
    "18",	"Quebec __ Région de la Capitale-Nationale",	"Quebec __ Capitale-Nationale",
    "19",	"Quebec __ Région de la Mauricie et du Centre-du-Québec",	"Quebec __ Mauricie",
    "20",	"Quebec __ Région de l'Estrie",	"Quebec __ Estrie",
    "21",	"Quebec __ Région de Montréal",	"Quebec __ Montréal",
    "22",	"Quebec __ Région de l'Outaouais",	"Quebec __ Outaouais",
    "23",	"Quebec __ Région de l'Abitibi-Témiscamingue",	"Quebec __ Abitibi-Témiscamingue",
    "24",	"Quebec __ Région de la Côte-Nord",	"Quebec __ Côte-Nord",
    "25",	"Quebec __ Région du Nord-du-Québec",	"Quebec __ Nord-du-Québec",
    "26",	"Quebec __ Région de la Gaspésie - Îles-de-la-Madeleine",	"Quebec __ Gaspésie-Îles-de-la-Madeleine",
    "27",	"Quebec __ Région de la Chaudière-Appalaches",	"Quebec __ Chaudière-Appalaches",
    "28",	"Quebec __ Région de Laval",	"Quebec __ Laval",
    "29",	"Quebec __ Région de Lanaudière",	"Quebec __ Lanaudière",
    "30",	"Quebec __ Région des Laurentides",	"Quebec __ Laurentides",
    "31",	"Quebec __ Région de la Montérégie",	"Quebec __ Montérégie",
    "32",	"Quebec __ Région du Nunavik",	"Quebec __ Nunavik",
    "33",	"Quebec __ Région des Terres-Cries-de-la-Baie-James",	"Quebec __ Terres-Cries-de-la-Baie-James",
    "34",	"Ontario __ Erie St. Clair Health Integration Network", "",
    "35",	"Ontario __ South West Health Integration Network",	"",
    "36",	"Ontario __ Waterloo Wellington Health Integration Network",	"",
    "37",	"Ontario __ Hamilton Niagara Haldimand Brant Health Integration Network",	"",
    "38",	"Ontario __ Central West Health Integration Network",	"",
    "39",	"Ontario __ Mississauga Halton Health Integration Network",	"",
    "40",	"Ontario __ Toronto Central Health Integration Network",	"",
    "41",	"Ontario __ Central Health Integration Network",	"",
    "42",	"Ontario __ Central East Health Integration Network",	"",
    "43",	"Ontario __ South East Health Integration Network",	"",
    "44",	"Ontario __ Champlain Health Integration Network",	"",
    "45",	"Ontario __ North Simcoe Muskoka Health Integration Network",	"",
    "46",	"Ontario __ North East Health Integration Network",	"",
    "47",	"Ontario __ North West Health Integration Network",	"",
    "48",	"Ontario __ District of Algoma Health Unit",	"Ontario __ Algoma",
    "49",	"Ontario __ Brant County Health Unit",	"Ontario __ Brant",
    "50",	"Ontario __ Durham Regional Health Unit",	"Ontario __ Durham",
    "51",	"Ontario __ Grey Bruce Health Unit",	"Ontario __ Grey Bruce",
    "52",	"Ontario __ Haldimand-Norfolk Health Unit",	"Ontario __ Haldimand-Norfolk",
    "53",	"Ontario __ Halton Regional Health Unit",	"Ontario __ Halton",
    "54",	"Ontario __ City of Hamilton Health Unit",	"Ontario __ Hamilton",
    "55",	"Ontario __ Hastings and Prince Edward Counties Health Unit",	"Ontario __ Hastings Prince Edward",
    "56",	"Ontario __ Huron County Health Unit",	"",
    "57",	"Ontario __ Chatham-Kent Health Unit",	"Ontario __ Chatham-Kent",
    "58",	"Ontario __ Lambton Health Unit",	"Ontario __ Lambton",
    "59",	"Ontario __ Middlesex-London Health Unit",	"Ontario __ Middlesex-London",
    "60",	"Ontario __ Niagara Regional Area Health Unit",	"Ontario __ Niagara",
    "61",	"Ontario __ North Bay Parry Sound District Health Unit",	"Ontario __ North Bay Parry Sound",
    "62",	"Ontario __ Northwestern Health Unit",	"Ontario __ Northwestern",
    "63",	"Ontario __ City of Ottawa Health Unit",	"Ontario __ Ottawa",
    "64",	"Ontario __ Peel Regional Health Unit",	"Ontario __ Peel",
    "65",	"Ontario __ Perth District Health Unit",	"Ontario __ Huron Perth",
    "66",	"Ontario __ Peterborough County-City Health Unit",	"Ontario __ Peterborough",
    "67",	"Ontario __ Porcupine Health Unit",	"Ontario __ Porcupine",
    "68",	"Ontario __ Renfrew County and District Health Unit",	"Ontario __ Renfrew",
    "69",	"Ontario __ Eastern Ontario Health Unit",	"Ontario __ Eastern",
    "70",	"Ontario __ Simcoe Muskoka District Health Unit",	"Ontario __ Simcoe Muskoka",
    "71",	"Ontario __ Sudbury and District Health Unit",	"Ontario __ Sudbury",
    "72",	"Ontario __ Thunder Bay District Health Unit",	"Ontario __ Thunder Bay",
    "73",	"Ontario __ Timiskaming Health Unit",	"Ontario __ Timiskaming",
    "74",	"Ontario __ Waterloo Health Unit",	"Ontario __ Waterloo",
    "75",	"Ontario __ Wellington-Dufferin-Guelph Health Unit",	"Ontario __ Wellington Dufferin Guelph",
    "76",	"Ontario __ Windsor-Essex County Health Unit",	"Ontario __ Windsor-Essex",
    "77",	"Ontario __ York Regional Health Unit",	"Ontario __ York",
    "78",	"Ontario __ Oxford Elgin St. Thomas Health Unit",	"",
    "79",	"Ontario __ City of Toronto Health Unit",	"Ontario __ Toronto",
    "80",	"Manitoba __ Winnipeg Regional Health Authority",	"Manitoba __ Winnipeg",
    "81",	"Manitoba __ Prairie Mountain Health",	"Manitoba __ Prairie Mountain",
    "82",	"Manitoba __ Interlake-Eastern Regional Health Authority",	"Manitoba __ Interlake-Eastern",
    "83",	"Manitoba __ Northern Regional Health Authority",	"Manitoba __ Northern",
    "84",	"Manitoba __ Southern Health",	"Manitoba __ Southern Health",
    "85",	"Saskatchewan __ Sun Country Regional Health Authority",	"",
    "86",	"Saskatchewan __ Five Hills Regional Health Authority",	"",
    "87",	"Saskatchewan __ Cypress Regional Health Authority",	"",
    "88",	"Saskatchewan __ Regina Qu'Appelle Regional Health Authority",		"Saskatchewan __ Regina",
    "89",	"Saskatchewan __ Sunrise Regional Health Authority", "",
    "90",	"Saskatchewan __ Saskatoon Regional Health Authority",	"Saskatchewan __ Saskatoon",
    "91",	"Saskatchewan __ Heartland Regional Health Authority",	"",
    "92",	"Saskatchewan __ Kelsey Trail Regional Health Authority",	"",
    "93",	"Saskatchewan __ Prince Albert Parkland Regional Health Authority",	"",
    "94",	"Saskatchewan __ Prairie North Regional Health Authority",	"",
    "95",	"Saskatchewan __ Mamawetan/Keewatin/Athabasca",	"",
    "96",	"Saskatchewan __ Mamawetan Churchill River Regional Health Authority",	"",
    "97",	"Saskatchewan __ Keewatin Yatthé Regional Health Authority",	"",
    "98",	"Saskatchewan __ Athabasca Health Authority",	"",
    "99",	"Alberta __ South Zone",	"Alberta __ South",
    "100",	"Alberta __ Calgary Zone",	"Alberta __ Calgary",
    "101", "Alberta __ Central Zone",	"Alberta __ Central",
    "102",	"Alberta __ Edmonton Zone",	"Alberta __ Edmonton",
    "103",	"Alberta __ North Zone",	"Alberta __ North",
    "104",	"British Columbia __ East Kootenay Health Service Delivery Area",		"",
    "105",	"British Columbia __ Kootenay-Boundary Health Service Delivery Area",		"",
    "106",	"British Columbia __ Okanagan Health Service Delivery Area",		"",
    "107",	"British Columbia __ Thompson/Cariboo Health Service Delivery Area",		"",
    "108",	"British Columbia __ Fraser East Health Service Delivery Area",	"",
    "109",	"British Columbia __ Fraser North Health Service Delivery Area",	"",
    "110",	"British Columbia __ Fraser South Health Service Delivery Area",	"",
    "111",	"British Columbia __ Richmond Health Service Delivery Area",		"",
    "112",	"British Columbia __ Vancouver Health Service Delivery Area",	"BC __ Vancouver Coastal",
    "113",	"British Columbia __ North Shore/Coast Garibaldi Health Service Delivery Area",	"",
    "114",	"British Columbia __ South Vancouver Island Health Service Delivery Area",	"",
    "115",	"British Columbia __ Central Vancouver Island Health Service Delivery Area",	"",
    "116",	"British Columbia __ North Vancouver Island Health Service Delivery Area",	"",
    "117",	"British Columbia __ Northwest Health Service Delivery Area",	"",
    "118",	"British Columbia __ Northern Interior Health Service Delivery Area",		"",
    "119",	"British Columbia __ Northeast Health Service Delivery Area",	"",
    "120",		"",	"BC __ Not Reported",
    "121",		"",			"Alberta __ Not Reported",
    "122",		"",			"Nunavut __ Nunavut",
    "123",		"",			"NWT __ NWT",
    "124",		"",			"Yukon __ Yukon",
    "125",		"",		"PEI __ Prince Edward Island",
    "126",		"",		"Repatriated __ Not Reported",
    "127",		"",		"BC __ Interior",
    "128",		"",		"BC __ Northern",
    "129",		"",		"BC __ Island",
    "130",		"",		"Saskatchewan __ Central",
    "131",		"",		"Saskatchewan __ Far North",
    "132",		"",		"Saskatchewan __ South",
    "133",		"",		"Saskatchewan __ Not Reported",
    "134",		"",		"Saskatchewan __ North",
    "135",		"",		"Quebec __ Not Reported",
    "136",		"",		"Ontario __ Haliburton Kawartha Pineridge",
    "137",		"",		"Ontario __ Kingston Frontenac Lennox & Addington",
    "138",		"",		"Ontario __ Southwestern",
    "139",		"",		"BC __ Fraser",
    "140",		"",		"Ontario __ Leeds Grenville and Lanark"

    ), ncol=3, byrow=TRUE
  ))
  colnames( au_look ) = c("no", "statscan", "covid19")

  toreturn = NA

  if (returntable) toreturn = au_look

  if (!is.null(statscan)) {
    i = match( statscan, au_look$statscan )
    if (!is.na(i)) {
      toreturn = au_look$covid19[ i]
      if (toreturn =="" ) toreturn = NA
    }
  }
  if (!is.null(covid19)) {
    i = match( covid19, au_look$covid19 )
    if (!is.na(i)) {
      toreturn = au_look$statscan[i]
      if (toreturn =="" ) toreturn = NA
    }
  }


  return(toreturn)

}
