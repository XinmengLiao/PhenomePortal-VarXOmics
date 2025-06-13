exomiser='/mnt/storage_pool/Genomics/Exomiser/exomiser-cli-14.0.0'
input='/mnt/storage_pool/Genomics/Genome/test113/BC/'

java -jar $exomiser/exomiser-cli-14.0.0.jar \
  --analysis $input/sample.yml \
  --assembly hg38 \
  --exomiser.data-directory=$exomiser/data \
  --spring.config.location=$exomiser/application.properties
