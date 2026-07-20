source(file.path(Sys.getenv("SDY2583_REPO_ROOT", unset = "."), "R", "shared", "bootstrap.R"))

############################################################
# SDY2583 / ImmPort Blood Immunotypes
# CP24 Flow Cytometry Pilot Hazırlık ve İlk QC Scripti
#
# ARCHIVAL SCOPE NOTE:
# This is the verified CP24 pilot/codebook recovered from the project archive.
# It is not represented as the final full-cohort CP24 production pipeline.
#
# Amaç:
# - ImmPort SDY2583 metadata dosyalarını okumak
# - Subject-level yaş, grup ve cinsiyet bilgisini temizlemek
# - FCS manifestini düzenlemek
# - Panel bilgisini çıkarmak
# - CP24 panelinden dengeli pilot set hazırlamak
# - İndirilen FCS dosyalarını R'a tanıtmak
# - flowCore ile ilk FCS okuma, marker kontrolü, transformasyon ve QC yapmak
#
# Not:
# - Bu script bugüne kadar birlikte kullandığımız kodları toparlar.
# - Local paths are configured through config/paths.R or environment variables.
# - Her önemli adımın gerekçesi # ile açıklanmıştır.
############################################################


############################################################
# 0. Çalışma ortamını kaydetme / geri yükleme
############################################################

# RStudio'yu kapatmadan önce mevcut objeleri kaydetmek için:
# Bu komut o anki R ortamındaki objeleri .RData dosyasına yazar.
# Böylece RStudio kapansa bile ertesi gün kaldığın yerden devam edebilirsin.
# save.image(file = file.path(sd_codebook_dir(), "SDY2583_R_workspace.RData"))

# Daha önce kaydettiysen, yeniden yüklemek için:
# load(file.path(sd_codebook_dir(), "SDY2583_R_workspace.RData"))

# Ortamda hangi objeler var görmek için:
# ls()


############################################################
# 1. ImmPort indirme klasörünü ve Tab klasörünü tanıtma
############################################################

# Ana indirme klasörünü tanıtıyoruz.
# Burada path örnektir. Sende klasör adı farklıysa değiştir.
root_dir <- sd_immport_download_dir()
codebook_dir <- sd_codebook_dir()
dir.create(codebook_dir, recursive = TRUE, showWarnings = FALSE)

# Ana klasör var mı diye kontrol ediyoruz.
# TRUE çıkarsa yol doğru, FALSE çıkarsa klasör yolu yanlış demektir.
dir.exists(root_dir)

# root_dir içinde hangi dosyalar/klasörler var görelim.
list.files(root_dir)

# SDY2583 ana çalışma klasörüne giriyoruz.
base_dir <- file.path(root_dir, "SDY2583")

# base_dir gerçekten var mı?
dir.exists(base_dir)

# SDY2583 klasörü içinde hangi dosyalar var?
list.files(base_dir)

# Tab.zip açıldıktan sonra oluşan klasör genellikle "SDY2583-DR58_Tab" idi.
# Bunun içinde ayrıca "Tab" isimli bir alt klasör bulunuyordu.
tab_dir <- file.path(base_dir, "SDY2583-DR58_Tab", "Tab")

# Tab klasörü doğru mu?
dir.exists(tab_dir)

# Tab klasöründeki ilk dosya adlarını görelim.
list.files(tab_dir)[1:30]


############################################################
# 2. Subject-level yaş ve temel metadata okuma
############################################################

# arm_2_subject.txt dosyası subject-level metadata için kritik.
# Bu dosyada subject ID, yaş bilgisi ve arm bilgisi bulunuyor.
arm_subj <- read.delim(
  file.path(tab_dir, "arm_2_subject.txt"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Tablo boyutu: satır sayısı ve sütun sayısı.
dim(arm_subj)

# İlk birkaç satırı görüntüleme.
head(arm_subj)

# Sütun adlarını kontrol etme.
names(arm_subj)

# Toplam tekil subject sayısı.
length(unique(arm_subj$SUBJECT_ACCESSION))


############################################################
# 3. Yaş değişkenini temizleme
############################################################

# Yaş değişkeni MIN_SUBJECT_AGE sütununda kullanılabilir görünüyordu.
# 999 değeri gerçek yaş değil, missing/unknown kodu gibi duruyordu.
arm_subj$age <- arm_subj$MIN_SUBJECT_AGE

# 120 ve üzerindeki değerleri gerçekçi yaş kabul etmiyoruz.
# Özellikle 999 kodunu NA yapıyoruz.
arm_subj$age[arm_subj$age >= 120] <- NA

# Yaş dağılımını özetliyoruz.
summary(arm_subj$age)

# Kaç kişide yaş bilgisi var?
sum(!is.na(arm_subj$age))

# Kaç kişide yaş missing/999?
sum(is.na(arm_subj$age))

# Yaş aralığı.
range(arm_subj$age, na.rm = TRUE)


############################################################
# 4. Yaş gruplarını oluşturma
############################################################

# Pilot analiz için yaşları üç gruba ayırdık:
# Young_<40: 0-39
# Middle_40_59: 40-59
# Older_60plus: 60+
arm_subj$age_group <- cut(
  arm_subj$age,
  breaks = c(0, 39, 59, 120),
  labels = c("Young_<40", "Middle_40_59", "Older_60plus"),
  right = TRUE
)

# Yaş grubu dağılımı.
table(arm_subj$age_group, useNA = "ifany")


############################################################
# 5. Arm/cohort bilgisini okuma
############################################################

# arm_or_cohort.txt dosyasında arm kodlarının isimleri bulunuyor.
# Ancak bu dosya tek başına tüm cancer/healthy ayrımını vermedi.
arm <- read.delim(
  file.path(tab_dir, "arm_or_cohort.txt"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dim(arm)
head(arm)
names(arm)
print(arm)

# arm_subj ile arm bilgisini birleştiriyoruz.
arm_subj2 <- merge(
  arm_subj,
  arm,
  by = "ARM_ACCESSION",
  all.x = TRUE
)

dim(arm_subj2)
head(arm_subj2)
names(arm_subj2)

# SUBJECT_PHENOTYPE sütunu boştu.
# Bu nedenle cancer/healthy ayrımını buradan alamadık.
table(arm_subj2$SUBJECT_PHENOTYPE, useNA = "ifany")


############################################################
# 6. Hastalık/tanı ile ilişkili dosyaları arama
############################################################

# Tab klasöründe hastalık/tanı/phenotype ile ilgili dosyaların isimlerini arıyoruz.
# Not: list.files() içinde value = TRUE kullanılmaz; o grep() argümanıdır.
list.files(
  tab_dir,
  pattern = "disease|diagnosis|phenotype|clinical|medical|condition|cancer",
  ignore.case = TRUE
)

# Daha güvenli yöntem:
all_files <- list.files(tab_dir)

all_files[grepl(
  "disease|diagnosis|phenotype|clinical|medical|condition|cancer",
  all_files,
  ignore.case = TRUE
)]

# Lookup disease dosyasını okuyoruz.
lk_disease <- read.delim(
  file.path(tab_dir, "lk_disease.txt"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dim(lk_disease)
head(lk_disease, 30)
names(lk_disease)

# Çalışma düzeyi condition/disease dosyasını okuyoruz.
study_disease <- read.delim(
  file.path(tab_dir, "study_2_condition_or_disease.txt"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dim(study_disease)
head(study_disease)
names(study_disease)


############################################################
# 7. Tüm text dosyalarında cancer/healthy kelimelerini arama
############################################################

# Bu adımı cancer/healthy bilgisinin hangi dosyada geçtiğini bulmak için kullandık.
txt_files <- list.files(tab_dir, pattern = "\\.txt$", full.names = TRUE)

cancer_hits <- sapply(txt_files, function(f) {
  lines <- readLines(f, warn = FALSE)
  any(grepl(
    "cancer|carcinoma|melanoma|tumor|tumour|neoplasm|adenocarcinoma",
    lines,
    ignore.case = TRUE
  ))
})

basename(txt_files[cancer_hits])

healthy_hits <- sapply(txt_files, function(f) {
  lines <- readLines(f, warn = FALSE)
  any(grepl("healthy|control", lines, ignore.case = TRUE))
})

basename(txt_files[healthy_hits])


############################################################
# 8. Subject_2_Flow_cytometry_result dosyasını okuma
############################################################

# Cancer/healthy ayrımını asıl bu dosyada bulduk.
# Bu dosya ImmPort indirme ekranından ayrı olarak indirilmişti.
subject_flow_path <- list.files(
  root_dir,
  recursive = TRUE,
  pattern = "Subject_2_Flow|Flow_cytometry_result|cytometry_result",
  ignore.case = TRUE,
  full.names = TRUE
)

subject_flow_path

# Eğer birden fazla sonuç varsa .txt olanı seçiyoruz.
subject_flow_path <- subject_flow_path[grepl("\\.txt$", subject_flow_path)]

# Dosyayı okuyoruz.
subject_flow <- read.delim(
  subject_flow_path[1],
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dim(subject_flow)
head(subject_flow)
names(subject_flow)

# Bu dosyada boşluklu sütun adları vardı:
# "Subject Accession", "ARM Name", "File Name", "Original File Name" vb.


############################################################
# 9. Cancer patient / Healthy control dağılımını çıkarma
############################################################

# Flow/sample satırı düzeyinde grup dağılımı:
# 4420 cancer sample, 4080 healthy sample çıktı.
table(subject_flow$`ARM Name`)

# Tekil subject düzeyinde grup dağılımı:
# Cancer patient = 442, Healthy control = 408 çıktı.
tapply(
  subject_flow$`Subject Accession`,
  subject_flow$`ARM Name`,
  function(x) length(unique(x))
)


############################################################
# 10. Subject-level temiz demografi tablosunu oluşturma
############################################################

# subject_flow içinde her subject için birden fazla FCS satırı var.
# Önce subject_id + disease_group bilgisini tekilleştiriyoruz.
subject_group <- unique(subject_flow[, c("Subject Accession", "ARM Name")])

# arm_subj içindeki yaş bilgisiyle subject_group'u birleştiriyoruz.
subject_demo <- merge(
  arm_subj[, c("SUBJECT_ACCESSION", "age", "age_group")],
  subject_group,
  by.x = "SUBJECT_ACCESSION",
  by.y = "Subject Accession",
  all.x = TRUE
)

dim(subject_demo)
head(subject_demo)

# Grup dağılımı.
table(subject_demo$`ARM Name`, useNA = "ifany")

# Grup x yaş grubu dağılımı.
table(subject_demo$`ARM Name`, subject_demo$age_group, useNA = "ifany")

# Subject-level tabloda aynı subject birden fazla kez var mı?
anyDuplicated(subject_demo$SUBJECT_ACCESSION)

# Yaş özetini kontrol ediyoruz.
summary(subject_demo$age)


############################################################
# 11. Cinsiyet bilgisini subject-level tabloya ekleme
############################################################

# Cinsiyet subject_flow içinde "Gender" olarak vardı.
subject_sex <- unique(subject_flow[, c("Subject Accession", "Gender")])

# Cinsiyeti subject_demo'ya ekliyoruz.
subject_demo2 <- merge(
  subject_demo,
  subject_sex,
  by.x = "SUBJECT_ACCESSION",
  by.y = "Subject Accession",
  all.x = TRUE
)

# Cinsiyet dağılımı.
table(subject_demo2$Gender, useNA = "ifany")

# Grup x cinsiyet dağılımı.
table(subject_demo2$`ARM Name`, subject_demo2$Gender, useNA = "ifany")


############################################################
# 12. Panel bilgisini FCS dosya adından çıkarma
############################################################

# Original File Name içinde panel kodu vardı:
# DBG324_CP24.fcs gibi.
subject_flow$panel <- sub(".*_(CP[0-9]+)\\..*", "\\1", subject_flow$`Original File Name`)

# Panel dağılımı.
table(subject_flow$panel)

# Grup x panel dağılımı.
table(subject_flow$`ARM Name`, subject_flow$panel)

# Her subject için kaç panel var?
panel_per_subject <- tapply(
  subject_flow$panel,
  subject_flow$`Subject Accession`,
  function(x) length(unique(x))
)

summary(panel_per_subject)
table(panel_per_subject)


############################################################
# 13. Temiz minimal tabloları kaydetme
############################################################

# Subject-level minimal tablo.
subject_clean_minimal <- data.frame(
  subject_id = subject_demo2$SUBJECT_ACCESSION,
  age_years = subject_demo2$age,
  age_group = subject_demo2$age_group,
  disease_group = subject_demo2$`ARM Name`,
  sex = subject_demo2$Gender,
  stringsAsFactors = FALSE
)

# FCS manifest minimal tablo.
fcs_manifest_minimal <- data.frame(
  subject_id = subject_flow$`Subject Accession`,
  disease_group = subject_flow$`ARM Name`,
  sex = subject_flow$Gender,
  age_years = subject_flow$`Subject Age`,
  biosample_id = subject_flow$`Biosample Accession`,
  biosample_type = subject_flow$`Biosample Type`,
  experiment_id = subject_flow$`Experiment Accession`,
  expsample_id = subject_flow$`Expsample Accession`,
  result_schema = subject_flow$`Expsample Result Schema`,
  expsample_description = subject_flow$`Expsample Description`,
  reagent_name = subject_flow$`Reagent Name`,
  file_info_id = subject_flow$`File Info ID`,
  file_name = subject_flow$`File Name`,
  original_file_name = subject_flow$`Original File Name`,
  panel = subject_flow$panel,
  stringsAsFactors = FALSE
)

# CSV olarak kaydetme.
write.csv(
  subject_clean_minimal,
  file = file.path(codebook_dir, "SDY2583_subject_clean_minimal.csv"),
  row.names = FALSE
)

write.csv(
  fcs_manifest_minimal,
  file = file.path(codebook_dir, "SDY2583_fcs_manifest_minimal.csv"),
  row.names = FALSE
)


############################################################
# 14. Panel-marker haritası oluşturma
############################################################

# Her panel için reagent_name alanında marker + fluorochrome bilgileri vardı.
# Örn: CD8_BB700, CD3_BV605, CD57_BV786 gibi.
panel_reagents <- unique(subject_flow[, c("panel", "Reagent Name")])

parse_reagents <- function(panel_id, reagent_string) {
  x <- unlist(strsplit(reagent_string, ","))
  x <- trimws(x)
  x <- x[nchar(x) > 0]

  parts <- strsplit(x, "_")

  data.frame(
    panel = panel_id,
    marker_full = x,
    marker = sapply(parts, function(z) z[1]),
    fluorochrome = sapply(parts, function(z) paste(z[-1], collapse = "_")),
    stringsAsFactors = FALSE
  )
}

panel_marker_long <- do.call(
  rbind,
  lapply(seq_len(nrow(panel_reagents)), function(i) {
    parse_reagents(panel_reagents$panel[i], panel_reagents$`Reagent Name`[i])
  })
)

panel_marker_long <- unique(panel_marker_long)

marker_frequency <- aggregate(
  panel ~ marker,
  data = panel_marker_long,
  FUN = function(x) paste(sort(unique(x)), collapse = ", ")
)

marker_frequency$n_panels <- sapply(
  strsplit(marker_frequency$panel, ", "),
  length
)

panel_summary <- aggregate(
  marker ~ panel,
  data = panel_marker_long,
  FUN = function(x) paste(sort(unique(x)), collapse = ", ")
)

panel_summary$n_markers <- sapply(
  strsplit(panel_summary$marker, ", "),
  length
)

write.csv(panel_marker_long, file.path(codebook_dir, "SDY2583_panel_marker_long.csv"), row.names = FALSE)
write.csv(marker_frequency, file.path(codebook_dir, "SDY2583_marker_frequency.csv"), row.names = FALSE)
write.csv(panel_summary, file.path(codebook_dir, "SDY2583_panel_summary.csv"), row.names = FALSE)


############################################################
# 15. CP24 panelini seçme ve tüm CP24 manifestini oluşturma
############################################################

# CP24 panelini seçiyoruz.
# Bu paneli ilk pilot için seçtik çünkü T-cell differentiation/senescence eksenine uygundu.
cp24_all <- subject_flow[subject_flow$panel == "CP24", ]

# Subject-level yaş ve yaş grubu bilgisini ekliyoruz.
cp24_all <- merge(
  cp24_all,
  subject_demo[, c("SUBJECT_ACCESSION", "age", "age_group")],
  by.x = "Subject Accession",
  by.y = "SUBJECT_ACCESSION",
  all.x = TRUE
)

# CP24 toplam kayıt sayısı 850 olmalı.
nrow(cp24_all)

# Grup x yaş grubu dağılımı.
table(cp24_all$`ARM Name`, cp24_all$age_group, useNA = "ifany")

write.csv(
  cp24_all,
  file = file.path(codebook_dir, "SDY2583_CP24_all_subjects_manifest.csv"),
  row.names = FALSE
)


############################################################
# 16. Dengeli CP24 pilot set oluşturma
############################################################

# Pilot analizde her disease_group x age_group hücresinden 20 kişi seçtik.
# Böylece 2 grup x 3 yaş grubu x 20 kişi = 120 FCS dosyası oldu.
set.seed(42)

eligible <- cp24_all[
  !is.na(cp24_all$age_group) &
    cp24_all$`ARM Name` %in% c("Cancer patient", "Healthy control"),
]

groups <- split(
  eligible,
  list(eligible$`ARM Name`, eligible$age_group),
  drop = TRUE
)

cp24_pilot <- do.call(
  rbind,
  lapply(groups, function(df) {
    df[sample(seq_len(nrow(df)), size = min(20, nrow(df))), ]
  })
)

table(cp24_pilot$`ARM Name`, cp24_pilot$age_group)
nrow(cp24_pilot)

write.csv(
  cp24_pilot,
  file = file.path(codebook_dir, "SDY2583_CP24_balanced_pilot_120_manifest.csv"),
  row.names = FALSE
)

writeLines(
  cp24_pilot$`File Name`,
  con = file.path(codebook_dir, "SDY2583_CP24_pilot_file_names.txt")
)


############################################################
# 17. İndirilen FCS klasörünü R'a tanıtma
############################################################

# CP24 FCS folder is configured with SDY2583_CP24_FCS_DIR.
fcs_dir <- sd_fcs_dir("CP24")

# Klasör yolunu görelim.
fcs_dir

# Yol doğru mu?
dir.exists(fcs_dir)

# Klasör içindeki FCS dosyalarını listeleyelim.
downloaded_fcs <- list.files(
  fcs_dir,
  pattern = "\\.fcs$",
  ignore.case = TRUE,
  full.names = FALSE
)

length(downloaded_fcs)
head(downloaded_fcs)


############################################################
# 18. Pilot dosyalar eksiksiz mi kontrol etme
############################################################

pilot_file_names <- readLines(file.path(codebook_dir, "SDY2583_CP24_pilot_file_names.txt"))

missing_files <- setdiff(pilot_file_names, downloaded_fcs)

# 0 çıkarsa eksik dosya yok.
length(missing_files)
missing_files


############################################################
# 19. flowCore kurulumu ve yüklenmesi
############################################################

# flowCore Bioconductor paketidir.
# İlk kez kurarken aşağıdaki iki satır çalıştırılır.
# install.packages("BiocManager")
# BiocManager::install("flowCore")
#
# Eğer "Update all/some/none? [a/s/n]:" sorarsa n yazmak yeterlidir.
library(flowCore)


############################################################
# 20. İlk FCS dosyasını okuma
############################################################

pilot_fcs_paths <- file.path(fcs_dir, pilot_file_names)

length(pilot_fcs_paths)
all(file.exists(pilot_fcs_paths))

fcs1 <- read.FCS(
  pilot_fcs_paths[1],
  transformation = FALSE,
  truncate_max_range = FALSE
)

fcs1

# Event sayısı ve kanal sayısı.
dim(exprs(fcs1))

# Ham kanal adları.
colnames(exprs(fcs1))

# FCS parametre/marker haritası.
marker_map_cp24 <- pData(parameters(fcs1))
marker_map_cp24

write.csv(
  marker_map_cp24,
  file = file.path(codebook_dir, "SDY2583_CP24_marker_map_from_FCS.csv"),
  row.names = FALSE
)


############################################################
# 21. 120 pilot FCS için event count QC
############################################################

# Amaç:
# - 120 dosyanın tamamı okunuyor mu?
# - Her dosyada kaç event var?
# - Kanal sayısı hepsinde 19 mu?
event_qc <- data.frame(
  file_name = basename(pilot_fcs_paths),
  events = NA_integer_,
  n_channels = NA_integer_,
  read_ok = FALSE,
  error_message = NA_character_
)

for (i in seq_along(pilot_fcs_paths)) {
  cat("Reading", i, "of", length(pilot_fcs_paths), ":", basename(pilot_fcs_paths[i]), "\n")

  tmp <- tryCatch(
    read.FCS(
      pilot_fcs_paths[i],
      transformation = FALSE,
      truncate_max_range = FALSE
    ),
    error = function(e) e
  )

  if (inherits(tmp, "flowFrame")) {
    event_qc$events[i] <- nrow(exprs(tmp))
    event_qc$n_channels[i] <- ncol(exprs(tmp))
    event_qc$read_ok[i] <- TRUE
  } else {
    event_qc$error_message[i] <- tmp$message
  }
}

table(event_qc$read_ok)
summary(event_qc$events)
table(event_qc$n_channels)
head(event_qc)

write.csv(
  event_qc,
  file = file.path(codebook_dir, "SDY2583_CP24_pilot_event_QC.csv"),
  row.names = FALSE
)


############################################################
# 22. Compensation / spillover bilgisi var mı kontrol etme
############################################################

# FCS dosyasında compensation matrix var mı diye bakıyoruz.
keyword(fcs1)[grep("SPILL|SPILLOVER|COMP", names(keyword(fcs1)), ignore.case = TRUE)]


############################################################
# 23. İlk ham plotlar
############################################################

# Ham fluorescence değerleri çok geniş aralıkta olduğu için doğrudan plotlar okunmayabilir.
plot(
  exprs(fcs1)[, "FSC-A"],
  exprs(fcs1)[, "SSC-A"],
  pch = ".",
  xlab = "FSC-A",
  ylab = "SSC-A",
  main = "CP24 first file: FSC-A vs SSC-A"
)

plot(
  exprs(fcs1)[, "BV605-A"],
  exprs(fcs1)[, "PerCP-Cy5-5-A"],
  pch = ".",
  xlab = "CD3 / BV605-A",
  ylab = "CD8 / PerCP-Cy5-5-A",
  main = "CP24 first file: CD3 vs CD8"
)


############################################################
# 24. Logicle transformasyon
############################################################

# Flow cytometry fluorescence kanalları ham halde çok geniş aralıkta olur.
# Logicle transformasyon, negatif ve çok geniş pozitif değerleri daha okunabilir hale getirir.
param <- pData(parameters(fcs1))

fluoro_channels <- param$name[
  !grepl("FSC|SSC|Time", param$name)
]

# Kritik düzeltme:
# param$name seçildiğinde R, $P7N gibi isimleri koruyabiliyor.
# Bu isimler transform() sırasında hata yaratır.
# unname() ile vektörü isimsiz hale getiriyoruz.
fluoro_channels <- unname(fluoro_channels)

fluoro_channels

lgcl <- estimateLogicle(fcs1, channels = fluoro_channels)

fcs1_t <- transform(fcs1, lgcl)

dat_t <- as.data.frame(exprs(fcs1_t))

# Grafik için çok fazla event çizmek yerine 50.000 event örnekliyoruz.
set.seed(1)
idx <- sample(seq_len(nrow(dat_t)), size = min(50000, nrow(dat_t)))


############################################################
# 25. Logicle sonrası okunabilir plotlar
############################################################

plot(
  dat_t[idx, "BV605-A"],
  dat_t[idx, "PerCP-Cy5-5-A"],
  pch = ".",
  xlab = "CD3 / BV605-A, logicle",
  ylab = "CD8 / PerCP-Cy5-5-A, logicle",
  main = "CP24 first file: CD3 vs CD8, logicle transformed"
)

plot(
  dat_t[idx, "PerCP-Cy5-5-A"],
  dat_t[idx, "BV786-A"],
  pch = ".",
  xlab = "CD8 / PerCP-Cy5-5-A, logicle",
  ylab = "CD57 / BV786-A, logicle",
  main = "CP24 first file: CD8 vs CD57"
)

plot(
  dat_t[idx, "PerCP-Cy5-5-A"],
  dat_t[idx, "PE-CF594-A"],
  pch = ".",
  xlab = "CD8 / PerCP-Cy5-5-A, logicle",
  ylab = "PD1 / PE-CF594-A, logicle",
  main = "CP24 first file: CD8 vs PD1"
)

plot(
  dat_t[idx, "PerCP-Cy5-5-A"],
  dat_t[idx, "PE-A"],
  pch = ".",
  xlab = "CD8 / PerCP-Cy5-5-A, logicle",
  ylab = "CX3CR1 / PE-A, logicle",
  main = "CP24 first file: CD8 vs CX3CR1"
)

plot(
  dat_t[idx, "PE-Cy5-A"],
  dat_t[idx, "BB515-A"],
  pch = ".",
  xlab = "CD45RA / PE-Cy5-A, logicle",
  ylab = "CD27 / BB515-A, logicle",
  main = "CP24 first file: CD45RA vs CD27"
)

plot(
  dat_t[idx, "BB515-A"],
  dat_t[idx, "BV650-A"],
  pch = ".",
  xlab = "CD27 / BB515-A, logicle",
  ylab = "CD62L / BV650-A, logicle",
  main = "CP24 first file: CD27 vs CD62L"
)


############################################################
# 26. CP24 marker yorum notları
############################################################

# CP24 panelinde FCS marker map'inden gördüğümüz markerlar:
# CD3    = BV605-A
# CD8    = PerCP-Cy5-5-A
# CD62L  = BV650-A
# CD95   = BV711-A
# CD57   = BV786-A
# CD27   = BB515-A
# CX3CR1 = PE-A
# PD1    = PE-CF594-A
# CD45RA = PE-Cy5-A
# CXCR3  = BV421-A
# CXCR5  = PE-Cy7-A
# Dump/viability composite = BV510-A:
#   Viability_CD4_CD13_CD19_TCRgd
#
# Bu nedenle CP24 pilotunun ana ekseni:
# - CD8 T-cell differentiation
# - CD8 T-cell senescence-like phenotype
# - CD8 T-cell exhaustion/chronic activation
# - naive-like / central memory-like / effector-like / TEMRA-like ayrımlar
#
# İlk biyolojik readout adayları:
# - CD8+CD57+
# - CD8+PD1+
# - CD8+CX3CR1+
# - CD8+CD27-
# - CD8+CD62L-
# - CD8+CD45RA+CD27+CD62L+ naive-like
# - CD8+CD45RA+CD27-CD62L- TEMRA/terminal-like
# - CD8+CD27-CD62L- terminal/effector-like
#
# Model fikri:
# outcome ~ age + disease_group + age:disease_group + sex
# veya
# outcome ~ age_group + disease_group + sex
############################################################


############################################################
# 27. Workspace'i tekrar kaydetme
############################################################

# Analizin bu aşamasındaki tüm objeleri kaydetmek için:
# save.image(file = file.path(sd_codebook_dir(), "SDY2583_R_workspace_after_CP24_QC.RData"))
