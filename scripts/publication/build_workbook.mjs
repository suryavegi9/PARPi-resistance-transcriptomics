import fs from 'node:fs/promises';
import {createRequire} from 'node:module';
import {pathToFileURL} from 'node:url';
const require=createRequire(import.meta.url);
const {Workbook,SpreadsheetFile}=await import('@oai/artifact-tool');
const out=process.argv[2];
if(!out)throw new Error('Provide the isolated table assembly output directory.');
try { await fs.access(`${out}/Supplementary_Tables_S1-S3.xlsx`); throw new Error('Refusing overwrite'); }
catch(e) { if(e.code!=='ENOENT')throw e; }
const data=JSON.parse(await fs.readFile(`${out}/tables.json`,'utf8'));
const wb=Workbook.create();
function col(n){let s='';for(n++;n;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s;return s;}
for(const [name,table] of Object.entries(data)){
 const sh=wb.worksheets.add(name);sh.showGridLines=false;
 const last=col(table.columns.length-1), end=table.rows.length+1;
 const range=sh.getRange(`A1:${last}${end}`);range.values=[table.columns,...table.rows];
 range.format.font={name:'Arial',size:10,color:'#202020'};
 range.format.rowHeight=28;range.format.verticalAlignment='center';range.format.columnWidth=17;
 sh.getRange(`A1:${last}1`).format.fill='#DCE6F1';sh.getRange(`A1:${last}1`).format.font.bold=true;
 sh.getRange(`A1:${last}1`).format.wrapText=true;sh.getRange(`A1:${last}1`).format.rowHeight=65;
 sh.getRange(`A1:A${end}`).format.columnWidth=55;
 sh.getRange(`B1:B${end}`).format.columnWidth=28;
 for(let c=0;c<table.columns.length;c++){
  const field=table.columns[c];const r=sh.getRange(`${col(c)}2:${col(c)}${end}`);
  if(/padj/.test(field))r.setNumberFormat('0.000E+00');
  else if(/NES|delta|mean/.test(field)&&!(/significant|same_direction/.test(field)))r.setNumberFormat('0.000000');
  else if(/^n_|contributing/.test(field))r.setNumberFormat('0');
  if(field.includes('classification')&&c>1)sh.getRange(`${col(c)}1:${col(c)}${end}`).format.columnWidth=29;
 }
 for(let r=3;r<=end;r+=2)sh.getRange(`A${r}:${last}${r}`).format.fill='#F4F6F8';
 sh.freezePanes.freezeRows(1);sh.freezePanes.freezeColumns(2);
 const t=sh.tables.add(`A1:${last}${end}`,true,`${name}_data`);t.style='TableStyleLight1';
 sh.getRange(`A1:${last}1`).format.horizontalAlignment='center';
 sh.tabColor='#385D7D';
}
const notes=JSON.parse(await fs.readFile(new URL('./table_notes.json',import.meta.url),'utf8'));
for(const [name,rows] of Object.entries(notes)){
 const sh=wb.worksheets.add(name);sh.showGridLines=false;
 sh.getRange(`A1:C${rows.length}`).values=rows;
 const r=sh.getRange(`A1:C${rows.length}`);r.format.font={name:'Arial',size:10,color:'#202020'};r.format.wrapText=true;r.format.verticalAlignment='top';
 sh.getRange(`A1:A${rows.length}`).format.columnWidth=28;sh.getRange(`B1:B${rows.length}`).format.columnWidth=65;sh.getRange(`C1:C${rows.length}`).format.columnWidth=65;
 r.format.rowHeight=60;sh.getRange('A1:C1').format.rowHeight=26;sh.getRange('A1:C1').format.fill='#DCE6F1';sh.getRange('A1:C1').format.font.bold=true;
 for(let i=1;i<rows.length;i++){let max=Math.max(...rows[i].map((x,j)=>Math.ceil(String(x).length/[28,65,65][j])));sh.getRange(`A${i+1}:C${i+1}`).format.rowHeight=Math.max(32,max*15+12);}
 sh.freezePanes.freezeRows(1);
}
await wb.recalculate();
console.log((await wb.inspect({kind:'region',sheetId:'Table_S3',range:'A1:J4',maxChars:2200,tableMaxCols:10})).ndjson);
await (await SpreadsheetFile.exportXlsx(wb)).save(`${out}/Supplementary_Tables_S1-S3.xlsx`);
for(const name of [...Object.keys(data),...Object.keys(notes)]){
 const preview=await wb.render({sheetName:name,range:name.startsWith('Table')?'A1:H8':'A1:C6',scale:1.4,format:'png'});
 await fs.writeFile(`${out}/${name}_preview.png`,new Uint8Array(await preview.arrayBuffer()));
}
console.log('Workbook exported and previews rendered.');
