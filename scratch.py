import os
import json
import re
import pymupdf
import base64
import tiktoken
import traceback
import difflib 
from openai import OpenAI

os.chdir('projects/gpt-sc-review')

API_KEY = os.getenv('OPENAI_API_KEY')
CELL_TYPE = 'microglia'
MODEL_A = 'gpt-4.1-2025-04-14'
MODEL_B = 'gpt-4.1-mini-2025-04-14'
MAX_CTX_A = 1_047_576
MAX_CTX_B = 1_047_576
RESERVE_TOKENS_A = 32_768
RESERVE_TOKENS_B = 32_768
ENC = tiktoken.get_encoding('cl100k_base')

SYSTEM_PROMPT_A = '''
=============================================================
PROMPT A (High-Fidelity Summary)
=============================================================
SYSTEM MESSAGE:
You are an expert research assistant specialising in computational biology
and neuroscience. Your primary role is to summarise and analyse
**single-cell/nucleus RNA sequencing (sc/snRNA-seq studies)** in
neurological disorders, focusing on **one paper at a time**. You will
produce structured summaries that:

**Derive all scientific content strictly and exclusively from the paper**.

**Apply the tag system** (detailed below) to highlight priority findings,
confidence levels, and potential contradictions.

Use the instructions below to ensure **comprehensive coverage** of each
paper's **cell-type heterogeneity**, **disease mechanisms**, **validation
methods**, and **possible conflicts**—all while maintaining a **cohesive
narrative style**. Use general knowledge only to identify contradictions;
do **not** introduce new claims.

-------------------------------------------------------------
1) USER-SPECIFIED CELL TYPE FOCUS
-------------------------------------------------------------
- Focus on the cell type named by the variable CELL_TYPE in your
  application.
- Centre your summary around **how the chosen cell type is characterised**
  in the paper:
  - Systematically identify and describe **all distinct subtypes or cell
    states** of the `CELL_TYPE` reported in the paper. For **each**
    subtype, **explicitly list**:
    - Its specific designation or name used in the paper (e.g., Mic.1,
      Astro.A, DAM).
    - Its key **defining marker genes**, noting direction (up/down) if
      specified.
    - Its primary characteristics or proposed functional role (e.g.,
      homeostatic, inflammatory, stress-response, proliferative).
    - Its association with disease status, pathology (e.g., amyloid, tau),
      clinical variables, or experimental conditions.
  - Mention morphological or spatial findings relevant to that cell type if
    reported.
  - Include how these subtypes might fit into disease progression or aging
    trajectories.
  - Briefly note any **homeostatic** or baseline subpopulations if the
    paper differentiates them from disease-associated states.

-------------------------------------------------------------
2) LEVEL OF DETAIL & STYLE
-------------------------------------------------------------
- **Prioritise Biological Findings**: Keep paper identification, methods,
  and technicalities concise. **Expand** on:
  - Disease-associated vs. homeostatic subtypes,
  - Marker genes (with up/down regulation if stated),
  - Functional implications (lipid metabolism, inflammatory response,
    etc.),
  - Spatial or morphological validation (e.g., immunostaining, in situ
    hybridisation),
  - Aging/disease stage transitions (pseudotime, computational
    modelling),
  - Contradictions with other known models or prior data.
  - Host or genetic drivers (age, sex, APOE, GWAS variants) and
    quantitative activation/morphology scores (e.g., PAM stages,
    compactness) that influence the chosen cell type.
- **Soften causal claims** if based on cross-sectional or computational
  data (e.g., 'strongly associated with,' 'may contribute to,'
  'potentially mediates'). Use **cautious language** for causal or
  temporal claims unless there is direct experimental or longitudinal
  evidence.

-------------------------------------------------------------
3) ANALYSIS FRAMEWORK (XML-LIKE TAGS)
-------------------------------------------------------------
Use the following tags, but keep <metadata> and <methods> succinct.
**Focus** on <findings> and <clinical> details.

PAPER IDENTIFICATION <metadata>
- Full citation (Authors, Year, Journal)
- Disease focus
</metadata>

METHODOLOGY EXTRACTION <methods>
- Briefly note scRNA-seq vs snRNA-seq
- Tissue/region sampled
- Major steps or validation methods only if critical
</methods>

RESULTS ANALYSIS <findings>
- **Cell Type Proportions**: Any quantitative changes for the chosen cell
  type or subtypes (stats if provided).
- **Differential Gene Expression**: Key genes, direction of change,
  magnitude if reported.
- **Pathway Enrichment**: Summaries of significantly altered pathways
  (lipid metabolism, complement, etc.).
- **Cell Subtype Identification & Characterization**: Present a structured
  breakdown of **each** subtype identified for the chosen `CELL_TYPE`.
  For every subtype, detail:
  - Its **name/label** from the paper.
  - A list of its **key defining marker genes**.
  - A summary of its **functional signature** or associated pathways.
  - Its classification (e.g., homeostatic, disease-associated,
    intermediate).
  - Any significant **changes in its proportion** or specific
    associations (e.g., with disease stage, genotype, pathology load).
  - Relevant morphological or spatial validation data, if provided.
  - *Ensure this information is presented clearly for **each distinct
    subtype*** (e.g., using separate descriptive sentences or short
    paragraphs per subtype within the main summary).
- **Modulators & Metrics**: host or genetic factors (age, sex, risk
  alleles), quantitative activation scores or morphology metrics that
  modulate the cell type or its sub-states.
- **Gene Regulatory Networks**: Any relevant transcription
  factors/regulators.
- **Cell-Cell Communication**: Important ligand-receptor pairs or
  cross-talk.
- **Spatial Analysis**: If given, highlight in situ or morphological
  findings that validate subpopulations.
- **Aging/Disease Trajectories**: Note temporal modelling or
  stage-specific shifts, referencing how subtypes may evolve.
- **Genetic or Multi-omic Integration**: If eQTLs or other data link
  subtypes to risk variants or phenotypes, include them.
</findings>

DISEASE RELEVANCE <clinical>
- Disease-specific roles of the chosen cell type.
- Mechanistic insights (e.g., how subtypes drive or mitigate pathology)
  but **use guarded language** if purely associative.
- Possible therapeutic or biomarker implications.
</clinical>

-------------------------------------------------------------
4) OUTPUT FORMAT
-------------------------------------------------------------
Since each summary focuses on **one paper**, provide:

1) **Quick Reference (≈50–100 words)**
    - One or two sentences capturing **the most critical** findings about
      the chosen cell type (major subtypes, marker genes, disease
      impact).
    - **Include one clause naming any key genetic, demographic or
      pathological driver of the highlighted cell subtype(s)** (e.g.,
      'Mic.13 is APOE ε4-enriched').
    - Minimal mention of methodology unless essential.

2) **Detailed Summary (≈800-1000 words)**
    - Expand on the content outlined in Section 3 (<metadata>, <methods>,
      <findings>, <clinical>) in **paragraphs**. Avoid bullet lists.
    - Highlight morphological/spatial or temporal data if presented.
    - **Use the Tag System** (defined in Section 5) below to flag
      priorities, confidence levels, and potential contradictions within
      these paragraphs.

3) **Research Implications (≈100–200 words)**
    - Summarise open questions or next steps for the chosen cell type.
    - Mention whether subtypes or marker genes align with known
      classification schemes (if the paper references them).
    - **Flag any conflicts** with prior models or known data if relevant.

-------------------------------------------------------------
5) TAG SYSTEM & USAGE
-------------------------------------------------------------
- <keyFinding priority='1-3'> ... </keyFinding>
  - **priority='1'** for major novel subtypes **or strong host/genetic
    modifiers** of disease-associated states,
  - **priority='2'** for supporting or secondary results,
  - **priority='3'** for technical or minor points.

- <confidenceLevel>high | medium | low</confidenceLevel>
  - Base your rating on evidence strength (sample size, morphological
    validation).
  - If a finding is derived from computational or cross-sectional data,
    consider 'medium' confidence unless validated in vivo.

- <contradictionFlag>none | details</contradictionFlag>
  - **'details'** if the paper's claims conflict with known models or
    other studies.
  - Provide a brief note explaining the contradiction if possible.
  - If no conflicts are found, set <contradictionFlag>none</contradictionFlag>
    for each major claim.

**Important**:
- Insert these tags **within the Detailed Summary** paragraphs to
  highlight key points.
- Even if no contradictions exist, use
  <contradictionFlag>none</contradictionFlag> to confirm.

=============================================================
END OF PROMPT A
=============================================================
'''

SYSTEM_PROMPT_B = '''
=============================================================
PROMPT B (Fact-Checking & Refinement - No Count)
=============================================================
SYSTEM MESSAGE:
You are a **High-Fidelity Verification Agent**. Your sole purpose is to
ensure the **factual accuracy** of the provided `Initial Summary` by
verifying every claim against the `Full PDF Text` (and associated images).
You make **only essential corrections** while **preserving the original
structure and style absolutely**. Your focus is information pertaining
to the user-specified `CELL_TYPE`.

-------------------------------------------------------------
1) INPUTS
-------------------------------------------------------------
- `Full PDF Text & Images`: The complete text extracted from the
  original scientific paper and associated relevant images. This is
  the ***only*** source of truth.
- `Initial Summary`: A structured summary of the paper, focusing on
  the `CELL_TYPE`, generated by a previous process.
- `CELL_TYPE`: The specific biological cell type that is the focus
  of the summary verification.

-------------------------------------------------------------
2) CORE INSTRUCTIONS
-------------------------------------------------------------
- **Strict Source Grounding**:
  - Systematically examine each factual statement within the
    `Initial Summary`, especially those related to the `CELL_TYPE`.
  - For each statement, locate the ***exact supporting evidence***
    within the `Full PDF Text` or Figures.
  - Any statement lacking direct, unambiguous support from the
    source material **must be corrected or removed**.
  - Do ***not*** use external knowledge or make inferences beyond what
    is explicitly stated in the source material.

- **Minimal Correction Principle**:
  - **Priority 1: Correction:** Modify statements ***only*** to align
    them perfectly with the source material. Correct inaccuracies
    in data, relationships, terminology, gene names, etc. If a
    statement is entirely unsubstantiated, **remove it**.
  - **Priority 2: Preservation:** Retain the original wording,
    phrasing, and sentence structure whenever possible. **Avoid
    rephrasing** if the original is factually correct according
    to the source. Edits must be surgical.
  - **Priority 3: Constrained Addition:** Add information ***only***
    if its absence makes an existing, related statement about the
    `CELL_TYPE` factually ***incomplete or misleading*** according
    to the source material. Insert the ***minimum necessary***
    information directly from the source to rectify this.

- **Absolute Format Integrity**:
  - The output ***must*** replicate the ***exact*** structure and
    formatting of the `Initial Summary`.
  - Preserve the presence, order, and content (unless corrected) of:
    - Sections (e.g., Quick Reference, Detailed Summary, Research
      Implications, or others if present in the input).
    - Paragraph breaks, spacing, and overall layout.
    - All original XML-like tags, including attributes:
      (e.g., `<code><metadata></code>`, `<code><methods></code>`,
      `<code><findings></code>`, `<code><clinical></code>`,
      `<code><keyFinding priority='...'></code>`,
      `<code><confidenceLevel>...</code>`, `<code><contradictionFlag>...</code>`, etc.)
  - Do not alter tag attributes unless the correction requires it
    based ***only*** on evidence in the source material.

-------------------------------------------------------------
3) OUTPUT SPECIFICATION
-------------------------------------------------------------
- Generate ***only*** the full text of the revised summary, adhering
  strictly to the **Minimal Correction Principle** and **Absolute
  Format Integrity** rules.
- Do ***not*** include any other text, explanations, or comments.

=============================================================
END OF PROMPT B
=============================================================
'''

client = OpenAI(api_key=API_KEY)

def extract_text(pdf_path: str) -> str:
    doc = pymupdf.open(pdf_path)
    pages = [page.get_text() or '' for page in doc]
    doc.close()
    return '\f'.join(pages) 

def extract_text_and_images(pdf_path: str) -> tuple[str, list[dict]]:
    doc = pymupdf.open(pdf_path)
    full_text = []
    images_data = []
    min_dimension = 100
    max_images = 1500

    for page_num, page in enumerate(doc):
        if len(images_data) >= max_images:
            break
        full_text.append(page.get_text(sort=True) or '')
        img_list = page.get_images(full=True)
        
        for _, img in enumerate(img_list):
            if len(images_data) >= max_images:
                break 
            xref = img[0]
            try:
                base_image = doc.extract_image(xref)
                img_width = base_image.get("width", 0)
                img_height = base_image.get("height", 0)
                
                if img_width < min_dimension or img_height < min_dimension:
                    continue

                image_bytes = base_image["image"]
                ext = base_image["ext"]
                image_b64 = base64.b64encode(image_bytes).decode('utf-8')
                images_data.append({
                    "page": page_num + 1,
                    "b64": image_b64,
                    "ext": ext
                })
            except Exception as e:
                print(f'warning: could not extract image xref {xref} on '
                      f'page {page_num+1}: {e}')
                continue
    doc.close()
    cleaned_text = re.sub(r'\n\s*\n', '\n\n', '\f'.join(full_text)).strip()
    return cleaned_text, images_data

def count_tokens(text: str) -> int:
    return len(ENC.encode(text))

def generate_summary(pdf_text: str, images_data: list[dict], 
                     focus_cell_type: str) -> str:

    message_content = [{"type": "text", "text": pdf_text}]
    for img_data in images_data:
        message_content.append({
            "type": "text",
            "text": f"\n--- Figure Content (Page {img_data['page']}) ---"
        })
        message_content.append({
            "type": "image_url",
            "image_url": {
                "url": f"data:image/{img_data['ext']};base64,{img_data['b64']}",
                "detail": "low"
            }
        })
    resp = client.chat.completions.create(
        model=MODEL_A,
        messages=[
            {'role': 'system', 'content': SYSTEM_PROMPT_A},
            {'role': 'user', 'content': f"CELL_TYPE: {focus_cell_type}"},
            {'role': 'user', 'content': message_content},
        ],
        max_tokens=RESERVE_TOKENS_A,
        temperature=0.1,
        top_p=1.0,
        frequency_penalty=0.0,
        presence_penalty=0.0,
    )
    return resp.choices[0].message.content

def refine_summary(pdf_text: str, images_data: list[dict], 
                   initial_summary: str, 
                   focus_cell_type: str) -> str:

    content_list_main = []
    content_list_main.append({"type": "text", 
                             "text": f"Full PDF Text:\n{pdf_text}"})
    
    for img_data in images_data:
        content_list_main.append({
            "type": "text",
            "text": f"\n--- Figure Content (Page {img_data['page']}) ---"
        })
        content_list_main.append({
            "type": "image_url",
            "image_url": {
                "url": f"data:image/{img_data['ext']};base64,{img_data['b64']}",
                "detail": "low"
            }
        })
        
    content_list_main.append({
        "type": "text",
        "text": f"\n\nInitial Summary to Fact-Check:\n{initial_summary}"
    })

    resp = client.chat.completions.create(
        model=MODEL_B,
        messages=[
            {'role': 'system', 'content': SYSTEM_PROMPT_B},
            {'role': 'user', 'content': f'CELL_TYPE: {focus_cell_type}'},
            {'role': 'user', 'content': content_list_main},
        ],
        max_tokens=RESERVE_TOKENS_B,
        temperature=0.0,
        top_p=1.0,
        frequency_penalty=0.0,
        presence_penalty=0.0,
    )
    return resp.choices[0].message.content.strip()

def calculate_corrections(initial: str, final: str) -> int:
    d = difflib.Differ()
    diff = list(d.compare(
        initial.splitlines(keepends=True), 
        final.splitlines(keepends=True)))
    correction_blocks = 0
    in_block = False
    for line in diff:
        if line.startswith('+ ') or line.startswith('- '):
            if not in_block:
                correction_blocks += 1
                in_block = True
        elif line.startswith('  '):
            in_block = False        
    return correction_blocks

def load_summaries(json_path: str) -> dict:
    if not os.path.exists(json_path):
        return {}
    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            return json.load(f)
    except json.JSONDecodeError:
        print(f'warning: could not decode {json_path}, starting fresh.')
        return {}

def save_summaries(json_path: str, summaries: dict):
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(summaries, f, indent=2, ensure_ascii=False)

def main():
    in_dir = 'pdfs'
    out_dir = 'out'
    os.makedirs(out_dir, exist_ok=True)
    json_path = os.path.join(out_dir, f'{CELL_TYPE}_summaries.json')
    summaries = load_summaries(json_path)

    for pdf_file in sorted(os.listdir(in_dir)):
        pid = os.path.splitext(pdf_file)[0]
        if pid in summaries and summaries[pid].get('fact-checked summary'):
            print(f'[{pid}] exists with fact-checked summary. skipping.')
            continue
        pdf_path = os.path.join(in_dir, pdf_file)

        try:
            print(f'[{pid}] extracting text and images...', end=' ')
            pdf_text, images_data = extract_text_and_images(pdf_path)
            
            initial_summary = None
            if pid in summaries and summaries[pid].get('initial summary'):
                print(f'found existing initial summary...', end=' ')
                initial_summary = summaries[pid]['initial summary']
            else:
                print(f'generating summary...', end=' ')
                initial_summary = generate_summary(
                    pdf_text, images_data, CELL_TYPE
                )
                if pid not in summaries: summaries[pid] = {}
                summaries[pid].update({
                    'pdf': pid,
                    'cell_type': CELL_TYPE,
                    'initial summary': initial_summary,
                })
                save_summaries(json_path, summaries)
            
            if initial_summary is None:
                print(f"error: No initial summary found or generated for {pid}")
                continue
            print(f'refining summary...', end=' ')
            fact_checked_summary = refine_summary(
                pdf_text, images_data, initial_summary, CELL_TYPE
            )
            corrections = calculate_corrections(
                initial_summary, fact_checked_summary
            )
            summaries[pid].update({
                'fact-checked summary': fact_checked_summary,
                'corrections': corrections 
            })
            save_summaries(json_path, summaries) 
            print(f'done (Corrections: {corrections}).')

        except ValueError as e:
            print(f'skipped ({e})')
        except Exception as e:
            print(f'\nerror processing {pid}: {e}')
            traceback.print_exc()

def collect_summaries():
    out_dir = 'out'
    json_path = os.path.join(out_dir, f'{CELL_TYPE}_summaries.json')
    md_path = os.path.join(out_dir, f'{CELL_TYPE}_summaries.md')

    print(f'collecting summaries for {CELL_TYPE} into {md_path}...')
    summaries = load_summaries(json_path)
    
    if not summaries:
        print('no summaries found to collect.')
        return

    with open(md_path, 'w', encoding='utf-8') as md:
        for pid in sorted(summaries.keys()):
            data = summaries[pid]
            summary = data.get('fact-checked summary', '') 
            if summary:
                md.write(f'# summary for {pid} ({CELL_TYPE})\n\n')
                md.write(f'{summary}\n\n---\n\n')
            
    print('summary collection finished.')

if __name__ == '__main__':
    main()
    collect_summaries()