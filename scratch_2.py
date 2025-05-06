import os
import re
import json
import base64
import difflib
import traceback
import concurrent.futures
from collections import defaultdict
import pymupdf
import polars as pl
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
MAX_WORKERS = 3

PDF_INPUT_DIR = 'pdfs'
METADATA_PATH = 'metadata.csv'
OUTPUT_DIR = 'out'
JSON_OUTPUT_NAME = f'{CELL_TYPE}_summaries.json'
MD_OUTPUT_DIR_NAME = 'summaries_md'
COLLECTED_MD_FILENAME = f'{CELL_TYPE}_collected_summaries.md'
FILTERED_MD_FILENAME = f'{CELL_TYPE}_collected_summaries_filtered.md'

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
confidence levels, and potential contradictions / departures **explicitly
discussed within the paper**.

Use the instructions below to ensure **coverage** of each paper's
**cell-type heterogeneity**, **disease mechanisms**, **validation
methods**, and **possible conflicts**—all while maintaining a **cohesive
narrative style**. Do **not** introduce new claims based on general knowledge.

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
  - Contradictions or departures from prior data **explicitly discussed by
    the authors within the paper**.
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
**If the paper reports minimal or no significant findings for the `CELL_TYPE`
across most categories, summarize this overall pattern concisely in the
<findings> section introduction or relevant subsections, rather than
elaborating detailed negatives for every category.**

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

2) **Detailed Summary (≈800-1000 words, shorter if findings sparse)**
    - Expand on the content outlined in Section 3 (<metadata>, <methods>,
      <findings>, <clinical>) in **paragraphs**. Avoid bullet lists.
      **If findings for the `CELL_TYPE` are sparse, prioritize conciseness
      over the target word count.**
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
  - **'details'** **only if the paper *explicitly discusses* specific findings
    conflicting with known models or other studies.**
  - Provide a brief note explaining the comparison **as described by the
    authors.**
  - If no conflicts are found, set <contradictionFlag>none</contradictionFlag>
    for each major claim.

**Important**:
- Insert these tags **within the Detailed Summary** paragraphs to
  highlight key points. Ensure `<contradictionFlag>details</contradictionFlag>`
  reflects only explicitly discussed comparisons.
- Even if no contradictions exist, use
  <contradictionFlag>none</contradictionFlag> to confirm.

=============================================================
END OF PROMPT A
=============================================================
'''

SYSTEM_PROMPT_B = '''
=============================================================
PROMPT B (Fact-Checking & Structural Standardization)
=============================================================
SYSTEM MESSAGE:
You are a **High-Fidelity Verification and Standardization Agent**. Your
purpose is twofold:
1.  Ensure the **factual accuracy** of the provided `Initial Summary` by
    verifying every claim against the `Full PDF Text` (and associated images).
2.  Ensure the final output strictly adheres to a **Standard Output
    Structure**.

You apply corrections minimally and preserve original phrasing whenever
possible, focusing on information pertaining to the user-specified
`CELL_TYPE`.

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

**PRIORITY 1: STRICT SOURCE GROUNDING & FACTUAL CORRECTION**
- **Verify Facts:** Systematically examine each factual statement within the
  `Initial Summary`, especially those related to the `CELL_TYPE`. Locate the
  ***exact supporting evidence*** within the `Full PDF Text` or Figures.
- **Apply Minimal Correction Principle:**
  - **Correction:** Modify statements ***only*** to align them perfectly
    with the source material. Correct inaccuracies in data, relationships,
    terminology, gene names, etc.
  - **Removal:** If a statement is entirely unsubstantiated by the source
    material, **remove it**.
  - **Preservation:** Retain the original wording and phrasing whenever
    possible *if the statement is factually correct*. Edits must be surgical.
  - **Constrained Addition:** Add information ***only*** if its absence makes
    an existing, related statement about the `CELL_TYPE` factually
    ***incomplete or misleading*** according to the source. Insert the
    ***minimum necessary*** information directly from the source.
- **No External Knowledge:** Do ***not*** use external knowledge or make
  inferences beyond what is explicitly stated in the source material.

**PRIORITY 2: STRUCTURAL STANDARDIZATION**
- **Define Standard Output Structure:** The final output *must* follow this
  exact structure and order:
    ```
    **Quick Reference** (...)
    ---
    **Detailed Summary** (...)
    [This section MUST contain the following tags/sections internally,
     in order: <metadata>...</metadata>, <methods>...</methods>,
     <findings>...</findings>, <clinical>...</clinical>]
    ---
    **Research Implications** (...)
    ```
- **Enforce Structure:** After verifying/correcting the content (Priority 1),
  examine the structure of the (potentially modified) `Initial Summary`.
  - **If** the structure **matches** the Standard Output Structure, preserve
    it exactly (including paragraph breaks, spacing, tags, etc., as per
    Priority 1).
  - **If** the structure **deviates** (e.g., missing sections, incorrect
    order, extraneous sections like preliminary `<findings>` blocks):
    - **Reorganize:** Rearrange the *verified* content blocks to fit the
      Standard Output Structure.
    - **Remove Extraneous Content:** Delete any sections or content blocks
      from the `Initial Summary` that are not part of the Standard Output
      Structure (e.g., remove preliminary `<findings>` or `<clinical>`
      blocks that appear *before* the `**Quick Reference**`).
    - **Preserve During Reorganization:** While reorganizing, preserve the
      verified phrasing, sentences, paragraph breaks within moved blocks,
      and all original XML-like tags (e.g., `<keyFinding>`,
      `<confidenceLevel>`, `<contradictionFlag>`) associated with their
      content as much as possible. Do not alter tag attributes unless
      required by factual correction (Priority 1).

-------------------------------------------------------------
3) OUTPUT SPECIFICATION
-------------------------------------------------------------
- Generate ***only*** the full text of the revised summary.
- The output MUST adhere strictly to the **Minimal Correction Principle** for
  content and the **Standard Output Structure** for format.
- Do ***not*** include any other text, explanations, or comments about the
  changes made.

=============================================================
END OF PROMPT B
=============================================================
'''

SYSTEM_PROMPT_C = '''
=============================================================
PROMPT C (Identify Insufficient Detail)
=============================================================
SYSTEM MESSAGE:
You are a detail-oriented analysis assistant. Your task is to review a
collection of paper summaries provided in Markdown format. Each summary
is associated with a unique identifier (PID), typically found in the
header (e.g., "# summary for [PID] (...)").

Your goal is to identify summaries that offer **minimal specific detail**
about the biological `CELL_TYPE` specified by the user. "Minimal specific
detail" means the summary mentions the `CELL_TYPE` but provides very few
or no specifics about:
- Distinct subtypes or states reported in the paper.
- Key defining marker genes for those subtypes/states.
- Specific functional roles attributed to the `CELL_TYPE` or its
  subtypes in the context of the study (e.g., disease association,
  pathway involvement).
- Quantitative changes or significant findings directly related to the
  `CELL_TYPE`.

Summaries that simply state the `CELL_TYPE` was present or analysed,
without elaborating on its characteristics or findings based on the paper,
should be flagged.

-------------------------------------------------------------
INPUT:
-------------------------------------------------------------
1.  `Markdown Content`: A string containing multiple paper summaries,
    each marked with a header like "# summary for [PID] (...)".
2.  `CELL_TYPE`: The specific cell type focus (e.g., microglia).

-------------------------------------------------------------
OUTPUT SPECIFICATION:
-------------------------------------------------------------
- Return **only** a plain text list of the unique identifiers (PIDs)
  for the summaries identified as insufficient.
- Each PID should be on a **new line**.
- Do **not** include headers, explanations, or any other text.

Example Output:
PID_A
PID_C
PID_F

=============================================================
END OF PROMPT C
=============================================================
'''

SYSTEM_PROMPT_D = '''
=============================================================
PROMPT D (Group Synthesis) 
=============================================================
SYSTEM MESSAGE:
You are an expert scientific synthesizer specializing in computational biology
and neuroscience. Your task is to generate a cohesive synthesis document
about a specific `CELL_TYPE` within a defined `CONTEXT` (e.g., disease,
cross-disorder analysis). Your synthesis must be based **exclusively** on the
provided `Input Summaries`, each representing findings from a single
scientific paper identified by a unique Persistent Identifier (PID).

-------------------------------------------------------------
1. Inputs
-------------------------------------------------------------
* `CONTEXT`: The specific theme for this synthesis (e.g., "Alzheimer's
    Disease").
* `CELL_TYPE`: The biological cell type focus (e.g., "microglia").
* `Input Summaries`: A collection of structured text blocks. Each block
    contains verified findings from one paper and is clearly associated
    with its source `PID` (e.g., marked "PID: [Identifier]"). **These
    summaries are the sole source of information.**

-------------------------------------------------------------
2. Core Principles
-------------------------------------------------------------

A. Evidence Fidelity and Attribution Mandate:
   * **Strict Source Reliance:** All synthesized content MUST derive
       *exclusively* from the provided `Input Summaries`. Do NOT introduce
       any external knowledge, assumptions, or interpretations.
   * **Maintain Source Certainty:** Accurately reflect the level of
       certainty indicated in the summaries (e.g., use 'suggests',
       'associated with' if the summaries used such language).
   * **Meticulous & Mandatory PID Attribution:** This is CRITICAL. Every
       specific factual statement, finding, reported marker gene, pathway,
       characteristic, subtype description, experimental result, or
       conclusion presented MUST be followed immediately by the PID(s) of
       the source summary/summaries in the format `(PID1)` or `(PID1, PID3,
       PID7)`. **Attribute every piece of specific information;
       unattributed claims render the output invalid.**

B. Input-Driven Thematic Synthesis with Comparative Analysis:
   * **Goal:** Identify and synthesize the most significant themes,
       findings, and relationships regarding the `CELL_TYPE` in the
       `CONTEXT` that emerge *from the collective information within the
       `Input Summaries`*.
   * **Focus on Prominence:** Structure your synthesis around the topics
       and findings most strongly represented or consistently discussed
       across the provided summaries.
   * **Consensus & Divergence:** Actively compare findings *between*
       different PIDs. Explicitly highlight areas of consensus (similar
       findings reported by multiple sources) and divergence (differing
        findings, findings unique to specific studies, or potential
        inconsistencies).
   * **Example Themes (Look for if prominent, do not force):** While
       synthesizing, pay attention to potentially prominent themes like cell
       subtype heterogeneity (names, markers, functions, associations),
       relevant biological pathways/processes (inflammation, metabolism),
       genetic/host factor influences (APOE, TREM2, age, sex), and the
       overall inferred role of the `CELL_TYPE`. The importance given to
       these themes must reflect their prominence in the *input summaries*.

C. Emergent Markdown Structure:
   * **Structure Origin:** The organization of the output document must
       emerge logically and naturally from the prominent themes identified
       during synthesis (Principle B).
   * **No Fixed Template:** Do NOT force the content into a predefined
       structure or include sections for topics not substantially covered in
       the `Input Summaries`. Avoid empty sections.
   * **Logical Flow:** Use markdown headers (`##`, `###`) to create a clear
       and logical flow that best communicates the synthesized knowledge
       derived *from the inputs*.

-------------------------------------------------------------
3. Output Specifications
-------------------------------------------------------------
* Generate **only** the full markdown text of the synthesized document.
* Start the document **exactly** with a level 1 markdown header:
    `# Synthesis for [CONTEXT] ([CELL_TYPE])`.
* Use subsequent markdown headers (`##`, `###`) logically as determined by
    the emergent structure (Principle C).
* Ensure **all** specific claims and findings are followed by their source
    PID(s) in parentheses, per the Mandate (Principle A).
* Do **not** include any other text, explanations, apologies, or comments
    about the synthesis process.

=============================================================
END OF PROMPT D 
=============================================================
'''

client = OpenAI(api_key=API_KEY)

def extract_text_and_images(pdf_path: str) -> tuple[str, list[dict]]:
    doc = pymupdf.open(pdf_path)
    full_text = []
    images_data = []
    min_dimension = 100
    max_images = 1500
    for page_num, page in enumerate(doc):
        if len(images_data) >= max_images: break
        full_text.append(page.get_text(sort=True) or '')
        img_list = page.get_images(full=True)
        for _, img in enumerate(img_list):
            if len(images_data) >= max_images: break
            xref = img[0]
            base_image = doc.extract_image(xref)
            if base_image.get('width', 0) < min_dimension or \
               base_image.get('height', 0) < min_dimension:
                continue
            image_bytes = base_image['image']
            image_b64 = base64.b64encode(image_bytes).decode('utf-8')
            images_data.append({
                'page': page_num + 1,
                'b64': image_b64,
                'ext': base_image['ext']
            })
    doc.close()
    cleaned_text = re.sub(r'\n\s*\n', '\n\n', '\f'.join(full_text)).strip()
    return cleaned_text, images_data

def generate_openai_completion(model: str, system_prompt: str,
                               user_content: any, max_tokens: int,
                               temperature: float = 0.1):
    return client.chat.completions.create(
        model=model,
        messages=[
            {'role': 'system', 'content': system_prompt},
            {'role': 'user', 'content': user_content},
        ],
        max_tokens=max_tokens,
        temperature=temperature,
        top_p=1.0,
        frequency_penalty=0.0,
        presence_penalty=0.0,
    )

def generate_summary(pdf_text: str, images_data: list[dict],
                     focus_cell_type: str) -> str:
    message_content = [{'type': 'text', 'text': pdf_text}]
    for img_data in images_data:
        message_content.extend([
            {'type': 'text',
             'text': f"\n--- Figure (Page {img_data['page']}) ---"},
            {'type': 'image_url', 'image_url': {
                'url': f"data:image/{img_data['ext']};base64,{img_data['b64']}",
                'detail': 'high'}}
        ])
    user_prompt = [
        {'type': 'text', 'text': f'CELL_TYPE: {focus_cell_type}'},
        *message_content
    ]
    resp = generate_openai_completion(
        MODEL_A, SYSTEM_PROMPT_A, user_prompt, RESERVE_TOKENS_A
    )
    return resp.choices[0].message.content

def refine_summary(pdf_text: str, images_data: list[dict],
                   initial_summary: str, focus_cell_type: str) -> str:
    content_list_main = [
        {'type': 'text', 'text': f'Full PDF Text:\n{pdf_text}'}
    ]
    for img_data in images_data:
         content_list_main.extend([
            {'type': 'text',
             'text': f"\n--- Figure (Page {img_data['page']}) ---"},
            {'type': 'image_url', 'image_url': {
                'url': f"data:image/{img_data['ext']};base64,{img_data['b64']}",
                'detail': 'high'}}
        ])
    content_list_main.append({
        'type': 'text',
        'text': f'\n\nInitial Summary to Fact-Check:\n{initial_summary}'
    })
    user_prompt = [
        {'type': 'text', 'text': f'CELL_TYPE: {focus_cell_type}'},
        *content_list_main
    ]
    resp = generate_openai_completion(
        MODEL_B, SYSTEM_PROMPT_B, user_prompt, RESERVE_TOKENS_B,
        temperature=0.0
    )
    return resp.choices[0].message.content.strip()

def calculate_similarity_score(initial: str, final: str) -> float:
    norm = lambda text: '\n'.join(
        line.strip() for line in text.splitlines() if line.strip()
    )
    return difflib.SequenceMatcher(None, norm(initial), norm(final)).ratio()

def load_json(json_path: str) -> dict:
    with open(json_path, 'r', encoding='utf-8') as f:
        return json.load(f)

def save_json(json_path: str, data: dict):
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

def process_single_paper(pid: str, pdf_path: str,
                         existing_initial_summary: str | None) -> dict | None:
    print(f'[{pid}] extracting...')
    pdf_text, images_data = extract_text_and_images(pdf_path)
    initial_summary = existing_initial_summary
    if initial_summary is None:
        print(f'[{pid}] generating summary...')
        initial_summary = generate_summary(pdf_text, images_data, CELL_TYPE)
    else:
         print(f'[{pid}] using existing initial summary.')
    print(f'[{pid}] refining summary...')
    fact_checked_summary = refine_summary(
        pdf_text, images_data, initial_summary, CELL_TYPE
    )
    similarity = calculate_similarity_score(
        initial_summary, fact_checked_summary
    )
    print(f'[{pid}] finished (similarity: {similarity:.4f}).')
    return {
        'pdf': pid, 'cell_type': CELL_TYPE,
        'initial summary': initial_summary,
        'fact-checked summary': fact_checked_summary,
        'similarity score': similarity
    }

def run_phase0(in_dir: str, out_dir: str, json_name: str) -> dict:
    os.makedirs(out_dir, exist_ok=True)
    json_path = os.path.join(out_dir, json_name)
    summaries = load_json(json_path) if os.path.exists(json_path) else {}
    tasks = []
    pdf_files = sorted(
        f for f in os.listdir(in_dir) if f.lower().endswith('.pdf')
    )
    for pdf_file in pdf_files:
        pid = os.path.splitext(pdf_file)[0]
        if pid in summaries and summaries[pid].get('fact-checked summary'):
            continue
        pdf_path = os.path.join(in_dir, pdf_file)
        existing = summaries.get(pid, {}).get('initial summary')
        tasks.append((pid, pdf_path, existing))
    if not tasks:
        print('Phase 0: no new papers.')
        return summaries
    print(f'Phase 0: processing {len(tasks)} papers '
          f'({MAX_WORKERS} workers)...')
    results = {}
    executor_opts = {'max_workers': MAX_WORKERS}
    with concurrent.futures.ThreadPoolExecutor(**executor_opts) as executor:
        future_map = {
            executor.submit(process_single_paper, pid, path, initial): pid
            for pid, path, initial in tasks
        }
        for future in concurrent.futures.as_completed(future_map):
            pid = future_map[future]
            result = future.result()
            if result: results[pid] = result
    if results:
        print(f'\nPhase 0: successfully processed {len(results)} papers.')
        summaries.update(results)
        print(f'Phase 0: saving summaries to {json_path}...')
        save_json(json_path, summaries)
    else:
        print('\nPhase 0: no papers successfully processed.')
    return summaries

def create_md_summaries(json_path: str, md_dir: str) -> bool:
    print(f"\nCreating individual .md summaries in {md_dir}...")
    summaries_data = load_json(json_path)
    os.makedirs(md_dir, exist_ok=True)
    count = 0
    for pid, data in summaries_data.items():
        summary = data.get('fact-checked summary', '')
        if summary:
            md_path = os.path.join(md_dir, f"{pid}.md")
            with open(md_path, 'w', encoding='utf-8') as f:
                cell_type = data.get('cell_type', CELL_TYPE)
                f.write(f"# Summary for {pid} ({cell_type})\n\n")
                f.write(f"**Quick Reference**\n[Placeholder]\n\n---\n\n")
                f.write(f"**Detailed Summary**\n{summary}\n\n---\n\n")
                f.write(f"**Research Implications**\n[Placeholder]\n")
            count += 1
    print(f"Created {count} .md files.")
    return count > 0

def extract_detailed_summary_from_md(md_content: str) -> str:
    match = re.search(
        r"\*\*Detailed Summary\*\*\s*\n(.*?)(?=\n---\s*\n|\Z)",
        md_content, re.DOTALL | re.IGNORECASE
    )
    if match: return match.group(1).strip()
    match_no_sep = re.search(
        r"\*\*Detailed Summary\*\*\s*\n(.*)", md_content,
        re.DOTALL | re.IGNORECASE
    )
    return match_no_sep.group(1).strip() if match_no_sep else ""

def generate_group_synthesis(group_name: str, pids: list[str],
                             cell_type: str, summary_dir: str,
                             out_dir: str, context_map: dict):
    print(f"\nPhase 1: Processing {group_name} ({len(pids)} summaries)")
    summary_blocks = []
    processed_pids = []
    for pid in pids:
        md_path = os.path.join(summary_dir, f"{pid}.md")
        if os.path.exists(md_path):
            with open(md_path, 'r', encoding='utf-8') as f: content = f.read()
            detail = extract_detailed_summary_from_md(content)
            if detail:
                summary_blocks.append(f"PID: {pid}\n\n{detail}")
                processed_pids.append(pid)
    input_text = "\n\n---\n\n".join(summary_blocks)
    context_desc = context_map.get(group_name, group_name)
    print(f"  Sending {len(processed_pids)} summaries for synthesis...")
    user_prompt = (
        f"CONTEXT: {context_desc}\nCELL_TYPE: {cell_type}\n\n"
        f"Input Summaries:\n{input_text}"
    )
    resp = generate_openai_completion(
        MODEL_A, SYSTEM_PROMPT_D, user_prompt, RESERVE_TOKENS_A
    )
    synthesis_content = resp.choices[0].message.content
    out_path = os.path.join(out_dir, f"{cell_type}_{group_name}_synthesis.md")
    with open(out_path, 'w', encoding='utf-8') as f: f.write(synthesis_content)
    print(f"  Synthesis saved to: {out_path}")

def run_phase1(md_dir: str, out_dir: str, meta_path: str):
    os.makedirs(out_dir, exist_ok=True)
    meta_df = pl.read_csv(meta_path)\
        .with_columns(
            (pl.col('first_author') + " " + pl.col('year').cast(pl.Utf8))
            .alias('PID'))
    pid_to_disease = dict(zip(meta_df['PID'], meta_df['disease']))
    group_defs = {
        'AD': "Alzheimer's Disease", 'PD': "Parkinson's Disease",
        'SZ': "Schizophrenia", 'MS': "Multiple Sclerosis",
    }
    context_map = {
        'AD': "Alzheimer's Disease", 'PD': "Parkinson's Disease",
        'SZ': "Schizophrenia", 'MS': "Multiple Sclerosis",
        'Cross-Disorder': "Cross-Disorder Studies"
    }
    synthesis_groups = ['AD', 'PD', 'SZ', 'MS', 'Cross-Disorder']
    grouped_pids = defaultdict(list)
    all_pids = {
        os.path.splitext(f)[0] for f in os.listdir(md_dir)
        if f.lower().endswith('.md')
    }
    print(f"\nPhase 1: Found {len(all_pids)} summary files in {md_dir}")
    for pid in all_pids:
        disease = pid_to_disease.get(pid)
        group = 'Minor Disease'
        if disease:
            exact_match = False
            for g_key, d_name in group_defs.items():
                if disease == d_name:
                    group = g_key
                    exact_match = True
                    break
            if not exact_match and 'Cross-Disorder' in disease:
                 group = 'Cross-Disorder'
        grouped_pids[group].append(pid)
    print("\nPhase 1: Summaries per Group:")
    total = 0
    for group, pids in grouped_pids.items():
        print(f"- {group}: {len(pids)}")
        total += len(pids)
    print(f"(Total: {total})")
    print("-" * 20)
    for group_name in synthesis_groups:
        if group_name in grouped_pids and grouped_pids[group_name]:
            generate_group_synthesis(
                group_name, grouped_pids[group_name], CELL_TYPE,
                md_dir, out_dir, context_map
            )
        else:
            print(f"\nPhase 1: Skipping {group_name}: No summaries.")
    print("\nPhase 1 processing complete.")
    minor_count = len(grouped_pids.get('Minor Disease', []))
    print(f"Minor Disease summaries ({minor_count}) not synthesized.")
    print(f"Synthesis outputs in: {out_dir}")

def collect_summaries_from_json(out_dir: str, cell_type: str, json_name: str,
                                output_filename: str,
                                pids_to_exclude: set[str] | None = None):
    json_path = os.path.join(out_dir, json_name)
    md_path = os.path.join(out_dir, output_filename)
    exclude_set = pids_to_exclude or set()
    print(f'\nCollecting summaries into {md_path}...')
    summaries_data = load_json(json_path)
    count = 0
    written = 0
    with open(md_path, 'w', encoding='utf-8') as md:
        for pid in sorted(summaries_data.keys()):
            summary = summaries_data[pid].get('fact-checked summary', '')
            if summary:
                count += 1
                if pid not in exclude_set:
                    md.write(
                        f'# summary for {pid} ({cell_type})\n\n'
                        f'{summary}\n\n---\n\n'
                    )
                    written += 1
    print(f'Collected {count} summaries from JSON, '
          f'wrote {written} to {md_path}.')
    if exclude_set: print(f'Excluded {len(exclude_set)} pids.')

def identify_insufficient_summaries(md_path: str, cell_type: str) -> list[str]:
    print(f'\nIdentifying insufficient summaries from {md_path}...')
    with open(md_path, 'r', encoding='utf-8') as f: md_content = f.read()
    user_prompt = f'CELL_TYPE: {cell_type}\n\nMarkdown Content:\n{md_content}'
    max_tokens = max(4096, RESERVE_TOKENS_B // 4)
    resp = generate_openai_completion(
        MODEL_B, SYSTEM_PROMPT_C, user_prompt, max_tokens, 0.0
    )
    pids = [
        p for p in resp.choices[0].message.content.strip().splitlines() if p
    ]
    print(f'Found {len(pids)} summaries potentially lacking detail.')
    return pids

if __name__ == '__main__':
    print("Starting Phase 0: Individual PDF Summary Generation...")
    json_full_path = os.path.join(OUTPUT_DIR, JSON_OUTPUT_NAME)
    summary_data = run_phase0(PDF_INPUT_DIR, OUTPUT_DIR, JSON_OUTPUT_NAME)
    print("\nPhase 0 Finished.")

    md_dir_full_path = os.path.join(OUTPUT_DIR, MD_OUTPUT_DIR_NAME)
    md_created = create_md_summaries(json_full_path, md_dir_full_path)

    if md_created:
        print("\nStarting Phase 1: Group Synthesis...")
        run_phase1(md_dir_full_path, OUTPUT_DIR, METADATA_PATH)
        print("\nPhase 1 Finished.")
    else:
        print("\nSkipping Phase 1: Issues creating .md files.")

    print("\nStarting Post-Analysis Steps...")
    if summary_data:
        collected_md_path = os.path.join(OUTPUT_DIR, COLLECTED_MD_FILENAME)
        collect_summaries_from_json(
            OUTPUT_DIR, CELL_TYPE, JSON_OUTPUT_NAME, COLLECTED_MD_FILENAME
        )
        insufficient = identify_insufficient_summaries(
            collected_md_path, CELL_TYPE
        )
        if insufficient:
            print("\nInsufficient PIDs:")
            for pid in insufficient: print(pid)
            filtered_md_path = os.path.join(OUTPUT_DIR, FILTERED_MD_FILENAME)
            collect_summaries_from_json(
                OUTPUT_DIR, CELL_TYPE, JSON_OUTPUT_NAME,
                FILTERED_MD_FILENAME, set(insufficient)
            )
        else:
            print("\nNo insufficient summaries identified.")
    else:
        print("\nNo summary data from Phase 0 for post-analysis.")

    print("\nFull script execution finished.")

meta = pl.read_csv(METADATA_PATH)
list(meta['disease'].unique())


