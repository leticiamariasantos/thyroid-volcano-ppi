# Cache versionado de APIs e fontes remotas

Cada fonte usa `data/cache/<fonte>/<versão>/`. Objetos binários grandes (`*.rds`) são
regeneráveis e não entram no Git; seus manifests, checksums e metadados textuais entram.
Uma etapa pode usar o cache somente após validar tipo, dimensões e checksum dos inputs.
Sem cache válido e sem resposta remota válida, o pipeline falha explicitamente.
