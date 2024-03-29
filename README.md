# UGW: Bouwsteen Oppervlaktewater

Python scripts voor oppervlaktewater-component Utrechts GrondWater model (UGW).

Voor een uitgebreide uitleg van de Bouwsteen Oppervlaktewater wordt verwezen 
naar het [rapport](./report/rap_bouwsteen_oppervlaktewater_UGW_incl_bijlage.pdf).

## Inhoudsopgave

-   [Installatie](#Installatie)
-   [Inhoud repository](#Inhoud-repository)
-   [Vragen, problemen of suggesties?](#Vragen-problemen-of-suggesties)
-   [Overzicht workflow](#Overzicht-workflow)
    -   [Ophalen oppervlaktewater gegevens](#Ophalen-oppervlaktewater-gegevens)
    -   [Omzetting naar MODFLOW-invoer](#Omzetting-naar-MODFLOW-invoer)
-   [Lokale data](#Lokale-data)
-   [Auteurs](#Auteurs)

## Installatie

De instructies hieronder gaan ervan uit dat `conda` (een Python package 
management tool) beschikbaar is. `conda` zit standaard in de Anaconda 
Python distributie.

1.  Als je nog geen Python op je computer hebt staan, 
    download [Anaconda](https://www.anaconda.com/distribution/).
2.  Clone of download de zip van deze repository en pak de bestanden uit.
3.  Open een terminal/command prompt en navigeer naar de map die je zojuist 
    uitgepakt hebt. Gebruik het volgende commando om een Python environment 
    te maken: `conda env create -f environment.yml`.
4.  Nu ben je klaar om de scripts te draaien! 
    Activeer de environment met `conda activate ugw`.

## Inhoud repository

De inhoud van de verschillende mappen in deze repository is hieronder kort 
toegelicht:

-   `config`: configuratie informatie voor toegang tot bepaalde webservices.
-   `data`: leeg, map voor het opslaan van alle benodigde data. Deze map wordt gevuld door 
    de verschillende scripts.
-   `models`: leeg, map voor het opslaan van modelinvoer
-   `scripts`: map met alle Python-scripts
    -   `oppervlaktewater`: map met scripts voor het ophalen en samenvoegen van oppervlaktewater gegevens
    -   `modflow`: map met scripts voor het omzetten van oppervlaktewater naar MODFLOW invoer
-   `tools`: Windows en Linux binaries voor MODFLOW
-   `environment.yml`: bestand met benodigde Python packages voor de scripts in deze repository.
-   `README.md`: dit bestand dat je nu leest

## Vragen, problemen of suggesties?

Als je een probleem tegenkomt in de code, of een suggestie hebt voor hoe het beter kan, of voor iets nieuws, dan kan dat via de [GitHub Issues pagina](https://github.com/ArtesiaWater/UGW/issues). Voordat je een issue meldt, kijk eerst of deze niet al door iemand anders is gemeld op deze pagina. Als dat niet zo is, dan kan je een [nieuwe issue](https://github.com/ArtesiaWater/UGW/issues/new/choose) aanmaken. Het development team zal op de issue reageren om het op te pakken. Bugs in de huidige code worden gefixt, voor nieuwe functionaliteit zal overlegd moeten worden.

Mocht je zelf de oplossing weten, of een bijdrage willen leveren aan de code, dan kan dat en dat wordt zeer op prijs gesteld. Zie [deze pagina ](https://github.com/ArtesiaWater/UGW/blob/8ded1cfbbbb017fa8f15ed0ab5a352ea705f1cbf/CONTRIBUTING.md) voor instructies hoe je suggesties of aanpassingen ter beoordeling kan klaarzetten door het core development team.

Vragen die niet direct over specifieke code/bugs/nieuwe features gaan kunnen gesteld worden via de [Discussions pagina](https://github.com/ArtesiaWater/UGW/discussions).

## Overzicht workflow

De scripts in deze repository voeren de volgende twee stappen uit:

1.  Downloaden en samenvoegen oppervlaktewatergegevens van verschillende 
    bronnen binnen het UGW interessegebied.
2.  Verwerken oppervlaktewater gegevens tot MODFLOW modelinvoer.

De volgorde waarin scripts moeten worden gedraaid is hieronder beschreven 
per stap. Daarin is ook aangegeven welke opties er beschikbaar zijn voor de 
gebruiker.

### Ophalen oppervlaktewater gegevens

Het ophalen, samenvoegen en valideren van oppervlaktewaterdata haalt de basisdata zoveel mogelijk uit online gegevensbronnen. Er zijn drie json-bestanden
aangemaakt die de gebruiker in principe niet hoeft aan te passen:

-   `config/administrations.json`: hier staan de organisaties (Rijkswaterstaat en waterschappen) waarvan data gedownloaded wordt, inclusief 
     een identificatie en een verwijzing naar een shape-file met een polygon per organisatie.
-   `config/sources.json`: hier staan per laag (watervlakken, waterlijnen en peilvakken) een verwijzing naar de web-service of het bestand waar de basisbestanden kan worden gedownload.
-   `validation.json`: hier staat per parameter een transformatie en defaultwaarde.

Helaas zijn niet alle data online beschikbaar en moet een deel lokaal beschikbaar gemaakt worden. Het bestand met deze data heet `sources_admin_waterschappen.zip` 
en kan [hier gedownload](https://drive.google.com/drive/u/0/folders/1fQ1pFgWc5xxuU5TBbR4kb75z1FNciWTj) worden. Het bestand dient uitgepakt te worden in de `data/sources` map.

Data ophalen voor een (deel)gebied kan door een sub-map te definieren in de data-map. In deze map
moet een shape-file staan met een of meerdere polygonen die samen het projectgebied definieren. De naam van deze map 
en de shape-file moeten worden opgegeven in het bestand config.ini in de map `scripts/oppervlaktewater`.

We gaan hieronder uit van de volgende structuur opgegeven in config.ini:

-   `data/project`: de map waarin de data wordt opgeslagen
-   `data/project/extent.shp`: polygon-shape met de grens van het project

De `scripts/config.ini` hoort er dan zo uit te zien:

    [general]
    project: project
    extend: extend.shp

De scripts zijn genummerd in volgorde:

-   `01_merge_sources.py`: download alle basisbestanden en voegt deze samen. Met bovenstaande `config.ini` in de map `data/project/input`
-   `02_prepare_modflow.py`: converteert en valideert de basisbestanden. Met bovenstaande `config.ini` in de map `data/project/modflow`

### Omzetting naar MODFLOW-invoer

De omzetting van de oppervlaktewater gegevens die in de vorige stap zijn 
aangemaakt naar MODFLOW invoer. De scripts zijn genummerd om aan te geven dat 
ze in een bepaalde volgorde uitgevoerd moeten worden.

-   `00_preprocessing.py`: zet benodigde data klaar in de data map
-   `01_create_model.py`: maak MODFLOW model met oppervlaktewater. De gebruiker 
    kan hier bepaalde dingen instellen
    -   `model_name` : modelnaam
    -   `extent` : gebied waarbinnen je een model wilt maken
    -   `delr` en `delc`: gridgrootte (dx en dy van het grid)
    -   aggregatiemethode voor oppervlaktewater (hoe moet oppervlaktewater omgezet worden naar modelinvoer). De opties zijn:
        -   `'individual'` : elke waterloop wordt individueel in het model opgenomen
        -   `'max_area'` : waterloop met grootste oppervlak binnen en cel wordt gebruikt voor bepaling van de parameters
        -   `'area_weighted'`: de parameters worden oppervlak-gewogen bepaald binnen en gridcel
        -   `'de_lange'`: de opschalingsformules van De Lange (1999) worden gebruikt voor de bepaling van de conductance van de waterlopen.
            Zie het rapport voor nadere toelichting.
-   `02*_postproc_*.py`: post-processing scrips, visualisatie van resultaten
-   `02_comparison.py`: vergelijking van modeluitkomsten tussen twee modellen.

## Lokale data

Niet alle data is online toegankelijk. Een deel moet lokaal gemaakt worden. 
Deze data kan via Google Drive gedownload worden. 
Zie daarvoor deze link: [Lokale data Google Drive](https://drive.google.com/drive/u/0/folders/1fQ1pFgWc5xxuU5TBbR4kb75z1FNciWTj)

Als je iets mist, kijk dan eerst of wat je zoekt via de link hierboven beschikbaar is.

## Auteurs

-   Daniel Tollenaar (D2Hydro)
-   Ruben Caljé (Artesia)
-   Davíd Brakenhoff (Artesia)
