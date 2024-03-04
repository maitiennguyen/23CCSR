import sys

argv = sys.argv

clustal_file = argv[1]
dataset_file = argv[2]

all_specs = """Saccharomyces paradoxus
Saccharomyces cerevisiae
Saccharomyces mikatae
Saccharomyces kudriavzevil
Saccharomyces uvarum
Saccharomyces eubayanus
Candida nivariensis
Nakaseomyces delphensis
Candida bracarensis
Candida glabrata
Nakaseomyces bacillisporus
Candida castellii
Kazachstania solicola
Kazachstania aerobia
Kazachstania unispora
Kazachstania siamensis
Kazachstania transvaalensis
Kazachstania yakushimaensis
Kazachstania taianensis
Kazachstania naganishii
Kazachstania bromeliacearum
Kazachstania intestinalis
Kazachstania martiniae
Kazachstania turicensis
Kazachstania kunashirensis
Kazachstania spencerorum
Kazachstania rosinii
Kazachstania africana
Kazachstania viticola
Naumovozyma dairenensis
Naumovozyma castelli
Tetrapisispora namnaonensis
Tetrapisispora fleeti
Tetrapisispora phaffii
Tetrapisispora iriomotensis
Vanderwaltozyma polyspora
Tetrapisispora blattae
Yueomyces sinensis
Torulaspora franciscae
Torulaspora pretoriensis
Torulaspora delbrueckli
Torulaspora maleeae
Torulaspora microellipsoides
Zygotorulaspora florentina
Zygotorulaspora mrakii
Zygosaccharomyces bisporus
Zygosaccharomyces bailii
Zygosaccharomyces kombuchaensis
Zygosaccharomyces rouxii
Lachancea lanzarotensis
Lachancea meyersii
Lachancea dasiensis
Lachancea nothofagi
Lachancea quebecensis
Lachancea thermotolerans
Lachancea waltii
Lachancea mirantina
Lachancea fermentati
Lachancea cidri
Lachancea kluyveri
Eremothecium gossypii
Ashbya aceri
Eremothecium cymbalariae
Eremothecium coryli
Eremothecium sinecaudum
Kluyveromyces lactis
Kluyveromyces dobzhanskii
Kluyveromyces marxianus
Kluyveromyces aestuarii
Hanseniaspora uvarum
Hanseniaspora nectarophila
Hanseniaspora clermontiae
Hanseniaspora meyeri
Hanseniaspora thailandica
Hanseniaspora opuntiae
Hanseniaspora pseudoguilliermondii
Hanseniaspora lachancei
Hanseniaspora guilliermondii
Kloeckera hatyaiensis
Hanseniaspora valbyensis
Hanseniaspora mollemarum
Hanseniaspora lindneri
Hanseniaspora jakobsenii
Hanseniaspora singularis
Hanseniaspora gamundiae
Hanseniaspora occidentalis
Hanseniaspora osmophila
Hanseniaspora vineae
Saccharomycodes ludwigii
Cyberlindnera jadinii
Candida vartiovaarae
Cyberlindnera maclurae
Cyberlindnera misumaiensis
Cyberlindnera suaveolens
Cyberlindnera mraki
Cyberlindnera americana
Candida mycetangii
Cyberlindnera petersonii
Cyberlindnera fabianii
Cyberlindnera xylosilytica
Candida freyschussii
Wickerhamomyces hampshirensis
Phaffomyces thermotolerans
Candida orba
Phaffomyces opuntiae
Phaffomyces antillensis
Candida montana
Barnettozyma pratensis
Barnettozyma salicaria
Barnettozyma californica
Barnettozyma populi
Barnettozyma hawaiiensis
Wickerhamomyces bovis
Wickerhamomyces canadensis
Wickerhamomyces sp.
Wickerhamomyces alni
Wickerhamomyces ciferri
Wickerhamomyces anomalus
Starmera amethionina
Candida stellimalicola
Starmera quercuum
Candida ponderosae
Ascoidea asiatica
Ascoidea rubescens
Saccharomycopsis crataegensis
Saccharomycopsis capsularis
Saccharomycopsis malanga
Metschnikowia matae var maris
Metschnikowia santaceciliae
Metschnikowia cerradonensis
Metschnikowia continentalis
Metschnikowia ipomoeae
Metschnikowia borealis
Metschnikowia hamakuensis
Metschnikowia mauinuiana
Metschnikowia hawaiiensis
Metschnikowia kamakouana
Metschnikowia dekortorum
Metschnikowia bowlesiae
Metschnikowia similis 
Metschnikowia arizonensis
Metschnikowia protea
Metschnikowia drakensbergensis
Metschnikowia shivogae
Metschnikowia aberdeeniae
Metschnikowia hibisci
Metschnikowia kipukae
Candida hawaiiana
Metschnikowia bicuspidata
Candida golubevii
Candida wancherniae
Candida intermedia
Candida blattae
Clavispora lusitaniae
Candida fructus
Candida oregonensis
Candida heveicola
Candida auris
Hyphopichia heimii
Candida rhagii
Candida gotoi
Danielozyma ontarioensis
Hyphopichia homilentoma
Hyphopichia burtonii
Candida parapsilosis
Candida orthopsilosis
Lodderomyces elongisporus
Candida corydali
Candida dubliniensis
Candida albicans
Candida tropicalis
Spathaspora arborariae
Spathaspora girioi
Spathaspora passalidarum
Spathaspora hagerdaliae
Spathaspora gorwiae
Scheffersomyces stipitis
Scheffersomyces lignosus
Suhomyces canberraensis
Suhomyces emberorum
Suhomyces tanzawaensis
Suhomyces pyralidae
Kodamaea ohmeri
Candida restingae
Kodamaea laetipori
Aciculoconidium aculeatum
Teunomyces cretensis
Teunomyces kruisii
Teunomyces gatunensis
Wickerhamia fluorescens
Priceomyces haplophilus
Priceomyces castillae
Priceomyces medius
Priceomyces carsonii
Debaryomyces subglobosus
Debaryomyces fabryi
Debaryomyces prosopidis
Debaryomyces nepalensis
Debaryomyces maramus
Millerozyma acaciae
Yamadazyma nakazawae
Yamadazyma philogaea
Yamadazyma tenuis
Yamadazyma scolyti
Candida gorgasii
Candida tammaniensis
Candida ascalaphidarum
Candida carpophila
Meyerozyma caribbica
Meyerozyma guilliermondii
Candida athensensis
Candida fragi
Kurtzmaniella cleridarum
Candida schatavii
Cephaloascus albidus
Cephaloascus fragrans
Babjeviella inositovora
Pichia norvegensis
Pichia kudriavzevii
Pichia occidentalis
Pichia exigua
Pichia heedii
Pichia nakasei
Pichia membranifaciens
Candida sorboxylosa
Pichia terricola
Saturnispora serradocipensis
Saturnispora hagleri 
Saturnispora dispora
Saturnispora zaruensis
Saturnispora saitoi
Saturnispora mendoncae
Saturnispora silvae
Kregervanrija fluxuum
Kregervanrija delftensis
Brettanomyces anomalus
Brettanomyces bruxellensis
Brettanomyces custersianus
Ambrosiozyma philentoma
Ambrosiozyma ambrosiae
Ambrosiozyma oregonensis
Ambrosiozyma monospora
Ambrosiozyma pseudovanderkliftii
Ambrosiozyma vanderkliftii
Ambrosiozyma kashinagacola
Ambrosiozyma maleeae
Ogataea polymorpha
Ogataea parapolymorpha
Ogataea philodendri
Ogataea kodamae
Ogataea minuta
Ogataea nonfermentans
Ogataea henricii
Ogataea pini
Ogataea glucozyma
Ogataea zsoltii
Ogataea trehaloabstinens
Ogataea trehalophila
Ogataea degrootiae
Ogataea methanolica
Candida arabinofermentans
Ogataea nitratoaversa
Ogataea pilisensis
Candida succiphila
Ogataea naganishii
Ogataea ramenticola
Candida boidinii
Kuraishia capsulata
Kuraishia molischiana
Kuraishia ogatae
Citeromyces hawaiiensis
Citeromyces matritensis
Citeromyces siamensis
Komagataella populi
Komagataella pseudopastoris
Komagataella mondaviorum
Komagataella pastoris
Komagataella ulmi
Komagataella phaffi
Komagataella kurtzmanii
Peterozyma xylosa
Peterozyma toletana
Nakazawaea peltata
Nakazawaea holstii
Pachysolen tannophilus 
Sporopachydermia lactativora 
Sporopachydermia quercuum 
Alloascoidea hylecoeti 
Blastobotrys raffinosifermentans 
Blastobotrys adeninivorans 
Blastobotrys peoriensis 
Blastobotrys americana 
Blastobotrys serpentis 
Blastobotrys proliferans 
Blastobotrys mokoenaii 
Blastobotrys muscicola 
Diddensiella caesifluorescens 
Spencermartinsiella europaea 
Sugiyamaella lignohabitans 
Zygoascus meyerae 
Zygoascus ofunaensis
Groenewaldozyma salmanticensis 
Starmerella apicola 
Starmerella bombicola 
Wickerhamiella domercqiae 
Wickerhamiella versatilis 
Wickerhamiella infanticola 
Wickerhamiella cacticola 
Deakozyma indianensis 
Candida incommunis
Dipodascus geniculatus
Dipodascus albidus
Geotrichum candidum 
Saprochaete clavata
Magnusiomyces tetraspermus 
Middelhovenomyces tepae 
Yarrowia divulgata 
Yarrowia deformans
Yarrowia lipolytica
Yarrowia keelungensis
Yarrowia bubula
Candida hispaniensis 
Nadsonia fulvescens var elongata 
Tortispora caseinolytica
Tortispora ganteri
Tortispora starmeri
Botryozyma nematodophila 
Trigonopsis vinaria
Trigonopsis variabilis
Lipomyces mesembrius 
Lipomyces arxii
Lipomyces tetrasporus
Lipomyces starkeyi
Lipomyces doorenjongii
Lipomyces japonicus 
Lipomyces lipofer
Lipomyces suomiensis
Lipomyces oligophaga""".splitlines()
all_specs = [spec.strip() for spec in all_specs]

blank_line_count = 0
rslt_specs = []
with open(clustal_file, 'r') as file:
    for line in file:
        if line.strip() == '':
            blank_line_count += 1
            if blank_line_count == 3:
                break
                
        item = line.split('\t')
        if len(item) > 1:
            spec = item[-1].strip()
            if spec not in rslt_specs:
                if "cteuk-1819" in spec:
                    spec = "Saccharomycetaceae sp."
                elif "cteuk-1144" in spec:
                    spec = "Saccharomycetales sp."
                elif "cteuk-1920" in spec:
                    spec = "Saccharomycopsidaceae sp."
                elif "Blastobotrys aristata" in spec:
                    spec = "Blastobotrys aristatus"
                elif "Saccharomyces kluyveri" in spec:
                    spec = "Lachancea kluyveri"
                elif spec == "Meyerozyma athensensis":
                    spec = "Candida athensensis"
                elif spec == "Kurtzmaniella fragi":
                    spec = 'Candida fragi'
                elif spec == "Clavispora fructus":
                    spec = 'Candida fructus'
                elif spec == "Hyphopichia gotoi":
                    spec = 'Candida gotoi'
                elif spec == "Metschnikowia hawaiiana": 
                    spec = 'Candida hawaiiana'
                elif spec == "Cyberlindnera mycetangii":
                    spec = 'Candida mycetangii'
                elif spec == "Hyphopichia rhagii":
                    spec = 'Candida rhagii'
                elif spec == "Starmera stellimalicola":
                    spec = 'Candida stellimalicola'
                elif spec == "Cyberlindnera mrakii":
                    spec = 'Cyberlindnera mraki'
                elif spec == "Grigorovia transvaalensis":
                    spec = 'Kazachstania transvaalensis'
                elif spec == "Grigorovia yakushimaensis":
                    spec = 'Kazachstania yakushimaensis'
                elif spec == "Hanseniaspora hatyaiensis":
                    spec = 'Kloeckera hatyaiensis'
                elif spec == "Kockiozyma suomiensis":
                    spec = 'Lipomyces suomiensis'
                elif spec == "Candida ipomoeae":
                    spec = 'Metschnikowia ipomoeae'
                elif spec == "Metschnikowia matae":
                    spec = 'Metschnikowia matae var maris'
                elif spec == "Metschnikowia proteae":
                    spec = 'Metschnikowia protea'
                elif spec == "Nadsonia fulvescens":
                    spec = "Nadsonia fulvescens var elongata"
                elif spec == 'Saccharomyces kudriavzevii':
                    spec = 'Saccharomyces kudriavzevil'
                elif spec == 'Tetrapisispora namnaoensis':
                    spec = 'Tetrapisispora namnaonensis'
                elif spec == 'Torulaspora delbrueckii':
                    spec = 'Torulaspora delbrueckli'
                elif spec == 'Wickerhamomyces ciferrii':
                    spec = 'Wickerhamomyces ciferri'
                elif spec == "Wickerhamomyces sp. nrrl yb-2243":
                    spec = "Wickerhamomyces sp."
                elif spec == 'Komagataella phaffii':
                    spec = 'Komagataella phaffi'
                elif spec == "Naumovozyma castellii":
                    spec = 'Naumovozyma castelli'
                elif spec == 'Tetrapisispora fleetii':
                    spec = 'Tetrapisispora fleeti'
                    
                if "sp." not in spec.lower() or "wickerhamomyces sp." == spec.lower():
                    rslt_specs.append(spec)


with open(dataset_file, 'a') as binary_file:
    for spec in all_specs:
        if spec in rslt_specs:
            binary_file.write('_'.join(spec.split()) + ",1\n")
        else:
            binary_file.write('_'.join(spec.split()) + ",0\n")
